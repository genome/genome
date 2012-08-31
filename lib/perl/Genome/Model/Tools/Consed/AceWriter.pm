package Genome::Model::Tools::Consed::AceWriter;

use strict;
use warnings;

use Tie::File;
use Data::Dumper;

class Genome::Model::Tools::Consed::AceWriter {
    has => [
        file => {
            is => 'Text',
            doc => 'Ace file to write',
        },
        width => {
            is => 'Integer',
            default_value => 50,
        },
        _fh => { is_optional => 1, },
        _contig_names => { default_value => {}, },
        _contig_count => { default_value => 0, },
        _read_count => { default_value => 0, },
        _at_buffer => { is_optional => 1, },
        _ct_buffer => { is_optional => 1, },
    ],
};

sub create {
    my $class = shift;

    my $self = $class->SUPER::create(@_);
    return if not $self;

    $self->_fh( IO::File->new('>'. $self->file) );
    $self->_fh->print("\n");
    $self->_at_buffer( IO::String->new() );
    $self->_ct_buffer( IO::String->new() );

    return $self;
}

sub close {
    my $self = shift;

    # contig tags
    $self->_ct_buffer->seek(0, 0);
    $self->_fh->print( join('', $self->_ct_buffer->getlines) );
    $self->_ct_buffer->close;

    # assembly tags
    $self->_at_buffer->seek(0, 0);
    $self->_fh->print( join('', $self->_at_buffer->getlines) );
    $self->_at_buffer->close;

    $self->_fh->close;
    
    # assembly header
    tie(my @array, 'Tie::File', $self->file);
    $array[0] = sprintf('AS %d %d', $self->_contig_count, $self->_read_count);
    untie(@array);

    return 1;
}

sub add_assembly_tags {
    my ($self, $tags) = @_;

    foreach my $tag ( @$tags ) {
        $self->_at_buffer->print (
            sprintf (
                "\nWA{\n%s %s%s%s%s\n}\n",,
                $tag->{tag_type},
                (exists $tag->{program}) ? $tag->{program}.' ' : '' ,
                $tag->{date},
                ( $tag->{no_trans} ? ' ' . $tag->no_trans : '' ),
                ( $tag->{data} ? "\n".$tag->{data} : '' ),
            )
        );
    }

    return 1;
}

sub add_contig_tags {
    my ($self, $tags) = @_;

    foreach my $tag ( @$tags ) {
        #print Dumper($tag);
        $self->_ct_buffer->print (
            sprintf (
                "CT{\n%s %s %s %d %d %s%s%s\n}\n\n",
                $tag->{contig_name},
                $tag->{tag_type},
                $tag->{program},
                $tag->{start_pos},
                $tag->{end_pos},
                $tag->{date},
                ( $tag->{no_trans} ? ' '.$tag->{no_trans} : ''),
                ( $tag->{data} ? "\n".$tag->{data} : '' ),
            )
        );
    }

    return 1;
}

sub add_contig { 
    my ($self, %p) = @_;

    my $contig = delete $p{contig};
    Carp::confess('No contig to export') if not $contig;
    $self->_increment_contig_count;

    # determine contig name
    my $name = ( $p{new_name} )
    ? delete $p{new_name}
    : ( $p{auto_rename} ) 
    ? sprintf('Contig%s', $self->_contig_count)
    : $contig->{name};
   
    # keep track of contig names
    Carp::confess("Contig with name ($name) already written to acefile") if exists $self->_contig_names->{$name};
    $self->_contig_names->{$name} = 1;

    my $read_count = $contig->{read_count};
    $self->_add_to_read_count($read_count);

    my $base_segments = $contig->{base_segments};
    my $base_segments_count = ( $contig->{base_segments} ) ? scalar (@$base_segments) : 0;

    $self->_fh->print (
        sprintf (
            "\nCO %s %d %d %d %s\n",
            $name,
            length $contig->{consensus},
            $contig->{read_count},
            $base_segments_count,
            $contig->{u_or_c},
        )
    );
    $self->_write_base_string($contig->{consensus});
    $self->_fh->print("\n");
    $self->_write_qualities($contig);
    $self->_write_read_positions($contig->{reads});
    $self->_write_base_segments($base_segments) if $contig->{base_segments};
    $self->_write_reads($contig->{reads});
    $self->_add_contig_tags_to_buffer($contig->{tags}) if $contig->{tags};
    
    return 1;
}

sub _increment_contig_count
{
    my $self = shift;

    return $self->_contig_count( $self->_contig_count + 1 );
}

sub _add_to_read_count
{
    my ($self, $count) = @_;

    return $self->_read_count( $self->_read_count + $count );
}

sub _write_base_string {
    my ($self, $base_string) = @_;
	
    my $width = $self->width;
	my $length = length($base_string);
    for (my $i = 0; $i < $length; $i += $width ) 
    {
        $self->_fh->print
        (
            substr
            (
                $base_string, 
                $i, 
                ( ( $i + $width - 1 ) < $length )
                ? $width 
                : $length - $i,
            ),
            "\n",
        );
    }

    return $self->_fh->print("\n");
}

sub _write_qualities 
{
	my ($self, $obj) = @_;

	$self->_fh->print("BQ\n");

    my $qualities = $obj->{base_qualities};
	my $width = $self->width;
	my $length = scalar(@$qualities) - 1;
    for (my $i = 0; $i <= $length; $i += $width ) {
        my $stop = ( $i + $width - 1 < $length )
        ? $i + $width - 1
        : $length;

        $self->_fh->print(' ', join(' ', @$qualities[$i..$stop]), "\n");
    }

    return $self->_fh->print("\n");
}

sub _write_read_positions {
    my ($self, $reads) = @_;

    for my $read ( sort { $a->{position} <=> $b->{position} || $a->{name} cmp $b->{name} } values %$reads ) {
        #print Dumper($read);
        $self->_fh->print (
            sprintf("AF %s %s %d\n", $read->{name}, $read->{u_or_c}, $read->{position}) 
        );
    }

    return 1;
}

sub _write_base_segments {
    my ($self, $base_segments) = @_;

    foreach my $bs ( @$base_segments ) {
        #print Dumper($bs);
        $self->_fh->print ( 
           sprintf("BS %d %d %s\n", $bs->{start_pos}, $bs->{end_pos}, $bs->{read_name})
        );
    }
    
    return $self->_fh->print("\n");
}

sub _write_reads {
    my ($self, $reads) = @_;

    for my $read ( sort { $a->{position} <=> $b->{position} || $a->{name} cmp $b->{name} } values %$reads ) {
        #print Dumper($read);
        my $tag_count = 0;
        my $tags = $read->{tags};
        $tag_count = ( $tags ) ? scalar @$tags : 0;
        $self->_fh->print (
            sprintf (
                "RD %s %d %d %d\n",
                $read->{name}, 
                length($read->{sequence}),
                $read->{info_count},
                $tag_count,
            )
        );
        
        $self->_write_base_string($read->{sequence});
        
        # QA
        $self->_fh->print (
            sprintf (
                "QA %d %d %d %d\n",
                $read->{qual_clip_start},
                $read->{qual_clip_end}, 
                $read->{align_clip_start},
                $read->{align_clip_end},
            )
        );

        # DS
        $self->_fh->print("DS");
        foreach my $attr (qw/ CHROMAT_FILE PHD_FILE TIME /) {
            my $val = $read->{description}->{$attr};
            next unless $val;
            $self->_fh->print (
                sprintf(' %s: %s', $attr, $val,)
            );
        }
        $self->_fh->print("\n");

        # Tags
        $self->_write_read_tags( $read->{tags} );
    }

    return 1;
}

sub _write_read_tags {
    my ($self, $tags) = @_;

    return unless $tags and @$tags;
    
    foreach my $tag ( @$tags ) {
        #print Dumper($tag);
        $self->_fh->print (
            sprintf
            (
                "RT{\n%s %s %s %d %d %s%s\n}\n\n",
                $tag->{read_name},
                $tag->{tag_type},
                $tag->{program},        
                $tag->{start_pos},
                $tag->{end_pos},
                $tag->{date},
                $self->{data} || '',
            )
        );
    }

    return 1;
}
   
1;

