package Genome::Model::Tools::Consed::AceToOut;

use strict;
use warnings;

use Genome;

use Bio::Seq::Quality;
use Bio::SeqIO;
use Data::Dumper;
use File::Basename;
use Genome::Sys; 

my @FORMATS = (qw/ fasta qual contig_names read_names /);
#fasta_for_tags => [qw/ get_bioseqs_for_tags fasta /],
#qual_for_tags => [qw/ get_bioseqs_for_tags qual /],
#tagstext => [qw/ get_tags tagstext /],
#bp_start_end => [qw/ get_tags bp_start_end /],
#bp_df => [qw/ get_tags  bp_tags /],
#bp_tags => [qw/ get_tags bp /],
#contignames => [qw/ contig_names contignames /],
my @BIOSEQ_FORMATS = (qw/ fasta qual /);

class Genome::Model::Tools::Consed::AceToOut {
    is => 'Command',
    has => [
    acefile => {
        is => 'Text',
        doc => 'Acefile',
    },
    format => {
        is => 'Text',
        doc => 'Format of output:' . join(', ', sort @FORMATS),
    },
    output_file => {
        is => 'Text',
        doc => 'Output file',
    },
    ],
    has_optional => [
    ctgs => {
        is => 'Text',
        doc => 'Specify contig names and poistions on the command line or in a file. Enter contig name separated by a comma (or one contig per line in a file) followed optionally by the unpadded positions.  Format the unpadded positions by following the contig name by an equals sign (=), then separating the start and stop positions by "to".  Example: Contig5 Contig50=33to1000 Contig77=1to100 This would get all of Contig5, Contig50 from 33 to 1000 and Contig77 from 1 to 100',
    },
    min_threshold => {
        is => 'Integer',
        doc => 'Contig minimum threshold, includes contigs with an unpadded length greater than or equal to this value',
    },
    max_threshold => {
        is => 'Integer',
        doc => 'Contig maximum threshold, includes contigs with an unpadded length less than or equal to this value',
    },
    inc_name => {
        is => 'Boolean',
        doc => 'Include acefile base name in output',
    },
    tag_types => {
        is => 'string',
        doc => 'If using a tag option, specify tag types to find, space delineated',
    },
    ],
};

sub help_brief {
    return '';
}

sub help_detail {
    return 'This script retrieves data from an acefile and then prints it directly to the screen.  To capture the output, use a redirect (>) to a file.';
}

sub formats {
    return @FORMATS;
}

sub execute {
    my $self = shift;

    $self->_validate_params
        or return;
    
    my $ace = Genome::Model::Tools::Consed::AceReader->create(file => $self->acefile)
        or return;

    my $filter_method = $self->{_filter_method};
    my $print_method = $self->{_print_method};
    while ( my $contig = $ace->next_contig ) {
        $self->$print_method($contig) if $self->$filter_method($contig);
    }

    $self->{_fh}->close;

    return 1;
}


sub Xexecute {
    my $self = shift;

    # Validate params
    $self->_validate_params
        or return;
    
    # Connect 
    my $ace = Genome::Model::Tools::Consed::AceReader->create(file => $self->acefile)
        or return;

    my $filter_method = $self->{_filter_method};
    my $print_method = $self->{_print_method};
    my $contig = $ace->next;
    CONTIG: while ( $contig ) {
        if ( $contig->{type} eq 'contig' ) {
            $contig->{unpadded_consensus} = $contig->{consensus};
            $contig->{unpadded_consensus} =~ s/\*//g;
            my $obj;
            my %read_positions;
            READ: while ( $obj = $ace->next ) {
                if ( $obj->{type} eq 'contig' ) {
                    last READ;
                }
                elsif ( $obj->{type} eq 'read_position' ) {
                    $read_positions{$obj->{read_name}} = $obj->{position}
                }
                elsif ( $obj->{type} eq 'read' ) {
                    $contig->{reads}->{ $obj->{name} } = $obj;
                    $contig->{reads}->{ $obj->{name} }->{position} = $read_positions{ $obj->{name} };
                }
            }
            $self->$print_method($contig) if $self->$filter_method($contig);
            $contig = $obj;
            next CONTIG;
        }
        $contig = $ace->next;
    }

    return $self->{_fh}->close;
}

sub _validate_params {
    my $self = shift;

    # Acefile
    Genome::Sys->validate_file_for_reading( $self->acefile )
        or return;

    # Format
    my $format = $self->format;
    unless ( grep { $format eq $_ } @FORMATS ) {
        $self->error_message("Invalid format ($format)");
        return;
    }
    $self->{_print_method} = '_print_' . $format;

    # Output
    $self->{_fh} = Genome::Sys->open_file_for_writing( $self->output_file )
        or return;
    $self->{_io} =  ( grep { $self->format eq $_ } @BIOSEQ_FORMATS ) 
    ? Bio::SeqIO->new('-format' => $format, '-fh' => $self->{_fh})
    : $self->{_fh};

    # Filtering
    my @threshold_properties = grep { defined $self->$_ } (qw/ max_threshold min_threshold /);
    if ( @threshold_properties and $self->ctgs ) { # Can't have threshold and ctgs
        $self->error_message("Can't filter by threshold and contig string");
        return;
    }
    # Threshold
    if ( @threshold_properties ) {
        for my $threshold_property ( @threshold_properties ) {
            next if $self->$threshold_property =~ /^\d+$/;
            $self->error_message("Threshold ($threshold_property) must be a positive integer");
            return;
        }
        if ( @threshold_properties == 1 ) {
            $self->{_filter_method} = '_filter_contig_by_' . $threshold_properties[0];
        }
        else {
            if ( $self->min_threshold > $self->max_threshold ) {
                $self->error_message("Invlaid thresholds.  The min is greater than the max");
                return;
            }
            $self->{_filter_method} = 'min_and_max_threshold';
        }
    }
    # Ctgs
    elsif ( my $ctgs = $self->ctgs ) {
        if ( -s $ctgs ) {
            my $fh = Genome::Sys->open_file_For_reading($ctgs)
                or return;
            $ctgs = join(' ', grep { /[\w\d]/ } grep { s/\s+//g } $fh->getlines);
            $fh->close;
        }
        my %ctgs;
        for my $ctg ( split(',', $ctgs) ) {
            my ($name, $pos) = split('=', $ctg);
            if ( $pos ) {
                $self->error_message("Invalid positions ($pos) in ctgs param for contig ($name)")
                    and return unless $pos =~ /^\d+to\d+$/;
                my ($start, $stop) = split('to', $pos);
                $self->{_ctgs}->{$name} = {
                    start => $start,
                    stop => $stop,
                };
            }
            else {
                $self->{_ctgs}->{$name} = undef;
            }

        }
        $self->{_filter_method} = '_filter_contig_by_ctgs';
    }
    else {
        $self->{_filter_method} = '_filter_contig_by_nothing';
    }

    return 1;
}

sub _filter_contig_by_nothing {
    return 1;
}

sub _filter_contig_by_max_and_min_threshold {
    my $self = shift;

    return ( $self->_filter_contig_by_max_threshold($_[0]) and $self->_filter_contig_by_min_threshold($_[0]) );
}

sub _filter_contig_by_max_threshold {
    my ($self, $contig) = @_;

    return length $contig->{unpadded_consensus} <= $self->max_threshold;
}

sub _filter_contig_by_min_threshold {
    my ($self, $contig) = @_;

    return length $contig->{unpadded_consensus} >= $self->min_threshold;
}

sub _filter_contig_by_ctgs {
    my ($self, $contig) = @_;

    return exists $self->{_ctgs}->{$contig->{name}};
}

sub _get_bioseq_from_contig {
    my ($self, $contig) = @_;

    my ($bases, $quals, $start, $stop);
    if ( my $pos = delete $self->{_ctgs}->{ $contig->{name} } ) {
        my $offset = $pos->{start} - 1;
        my $length = $pos->{stop} - $pos->{start} + 1;
        $bases = substr($contig->{unpadded_consensus}, $offset, $length);
        $quals = join(' ', splice(@{$contig->{base_qualities}}, $offset, $length));
    }
    else {
        $bases = $contig->{unpadded_consensus};
        $quals = join(' ', @{$contig->{base_qualities}});
        $start = 1;
        $stop = length $contig->{unpadded_consensus};
    }
    
    return Bio::Seq::Quality->new(
        '-id' => sprintf('%s%s', ( $self->inc_name ? $self->acefile : '' ), $contig->{name}),
        #'-desc' => "from $start to $stop",
        '-alphabet' => 'dna',
        '-seq' => $bases,
        '-qual' => $quals,
    );
}

sub _print_fasta {
    my $self = shift;
    
    return $self->{_io}->write_seq( $self->_get_bioseq_from_contig($_[0]) );
}

sub _print_qual {
    my $self = shift;

    return $self->{_io}->write_seq( $self->_get_bioseq_from_contig($_[0]) );
}

sub _print_contig_names {
    my ($self, $contig) = @_;
    
    return $self->{_io}->print( $contig->{name} . "\n" );
}

sub _print_read_names {
    my ($self, $contig) = @_;
    
    if ( my $pos = delete $self->{_ctgs}->{ $contig->{name} } ) {
        my @unpad_pos = $self->_get_unpad_positions_to_pad_positions_for_contig($contig);
        my $start = $unpad_pos[ $pos->{start} ];
        my $stop = $unpad_pos[ $pos->{stop} ];
        for my $read_name ( keys %{$contig->{reads}} ) {
            my $read = $contig->{reads}->{$read_name};
            next if $read->{position} > $stop;
            next if ($read->{position} + length($read->{sequence}) - 1) < $start;
            $self->{_io}->print( $read->{name} . "\n" );
        }
    }
    else{
        for my $read_name ( keys %{$contig->{reads}} ) {
            $self->{_io}->print( $read_name . "\n" );
        }
    }

    return 1;
}

sub _get_unpad_positions_to_pad_positions_for_contig {
    my ($self, $contig) = @_;

    my @pos;
    my $pad = 0;
    my $unpad = 0;
    foreach my $base ( split(//, $contig->{consensus}) ) {
        $pad++;
        $unpad++ if $base ne '*';
        $pos[$unpad] = $pad;
    }

    return @pos;
}

1;

