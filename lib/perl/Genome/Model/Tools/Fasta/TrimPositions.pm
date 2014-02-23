package Genome::Model::Tools::Fasta::TrimPositions;

use strict;
use warnings;

use Genome;
use Data::Dumper;

class Genome::Model::Tools::Fasta::TrimPositions {
    is  => 'Command',
    has => [
        fasta_file => {
            is => 'Text',
            doc => 'Input fasta file',
        },
        qual_file => {
            is => 'Text',
            doc => 'Input quality file',
            is_optional => 1,
        },
        trim_file => {
            is => 'Text',
            doc => 'File of trim positions, eg, Contig0.1 1-100,1000-end',
        },
        min_contig_length => {
            is => 'Number',
            doc => 'Exclude contigs less than this length',
            is_optional => 1,
            default => 1,
        },
        extend_trim_to_end_length => {
            is => 'Number',
            doc => 'Extend trimming to end if less than this many bases remaining',
            is_optional => 1,
            default => 0,
        },
    ],
};

sub help_brief {
    return 'Tool to trim fasta optionally with corrisponding qual at specified positions',
}

sub help_detail { 
    return <<EOS 
gmt fasta trim-position --fasta-file contigs.bases --qual-file contigs.quals --min-contig-length 200 --trim-file trim.txt

Example trim file:
 Contig0.1 start-100,2000-end
 Contig9.2 959-1000
 Contig9.2 1500-1550
EOS
}

sub execute {
    my $self = shift;

    $self->error_message('Input fasta is missing or empty: '.$self->fasta_file) and return if not
        -s $self->fasta_file;

    $self->error_message('Input trim-file is missing or empty: '.$self->trim_file) and return if not
        -e $self->trim_file;

    my $trim_pos;
    if ( not $trim_pos = $self->trim_positions_from_file ) {
        $self->error_message("Failed to get trim positions from file: ".$self->trim_file);
        return;
    }

    if ( not $self->trim_and_write_fasta($trim_pos) ) {
        $self->error_message("Failed to trim and write fasta file");
        return;
    }

    return 1;
}

sub trim_and_write_fasta {
    my ($self,$trim_pos) = @_;

    my %reader_params = (
        file => $self->fasta_file,
    );
    my %writer_params = (
        file => $self->fasta_file.'.clipped',
    );

    if ( $self->qual_file ) {
        $reader_params{qual_file} = $self->qual_file;
        $writer_params{qual_file} = $self->qual_file.'.clipped';
        unlink $self->qual_file.'.clipped';
    }

    my $reader = Genome::Model::Tools::Sx::PhredReader->create( %reader_params );
    unlink $self->fasta_file.'.clipped';
    my $writer = Genome::Model::Tools::Sx::PhredWriter->create( %writer_params );

    while ( my $seq = $reader->read ) {
        my $id = $seq->{id};
        my $seq_length = length $seq->{seq};
        if ( exists $trim_pos->{lc $id} ) {
            my @seqs = split( '', $seq->{seq} );
            my @quals;
            if ( $self->qual_file ) {
                @quals = split('', $seq->{qual});
            }
            for my $pos ( @{$trim_pos->{lc $id}} ) {
                my ( $start,$end ) = split ('-', $pos);
                # start trim position
                $start = ( $start =~ /start/i ) ? 1 : $start;
                $start = ( $start < 1 ) ? 1 : $start;
                $start = ( $start > $self->extend_trim_to_end_length ) ? $start : 1;
                # end trim position
                $end = ( $end =~ /end/i ) ? $seq_length : $end;
                $end = ( $end > $seq_length ) ? $seq_length : $end;
                $end = ( ($seq_length - $end) > $self->extend_trim_to_end_length ) ? $end : $seq_length;
                $self->debug_message("Trimming $id from $start to $end");
                # fasta
                for ( $start .. $end ) {
                    $seqs[$_ - 1] = '/';
                }
                # qual
                if ( $self->qual_file ) {
                    for ( $start .. $end ) {
                        undef $quals[$_ - 1];
                    }
                }
            }
            # remove clipped fasta
            my $clipped_seq = join('', @seqs);
            $clipped_seq =~ s/\///g;
            $seq->{seq} = $clipped_seq;
            # remove clipped qual
            if ( $self->qual_file ) {
                my $clipped_quals;
                for ( 0 .. $#quals ) {
                    next unless defined $quals[$_];
                    $clipped_quals .= $quals[$_];
                }
                $seq->{qual} = $clipped_quals;
            }
        }
        if ( length $seq->{seq} < $self->min_contig_length ) {
            $self->debug_message(
                "Removing $id (".(length $seq->{seq})." bp) is less than min length (".$self->min_contig_length.")"
            );
            next;
        }
        $writer->write( $seq );
    }

    return 1;
}

sub trim_positions_from_file {
    my $self = shift;

    my $fh = Genome::Sys->open_file_for_reading( $self->trim_file );
    my %pos;
    while ( my $line = $fh->getline ) {
        chomp $line;
        my ( $id ) = $line =~ /^(\S+)\s+/;
        my ( $p_string ) = $line =~ /^\S+\s+(.*)$/;
        $p_string =~ s/\s+//;
        my @ps = split(',', $p_string);
        for my $pos ( @ps ) {
            unless ( $pos =~ /^\d+-\d+$/ or $pos =~ /^start-\d+$/i or $pos =~ /^\d+-end/i or $pos =~ /^start-end$/i ) {
                $self->error_message("Expected position like 1-100 (START-END) but got: $pos");
                return;
            }
            push @{$pos{ lc $id }}, $pos;
        }
    }
    $fh->close;
    return \%pos;
}

1;
