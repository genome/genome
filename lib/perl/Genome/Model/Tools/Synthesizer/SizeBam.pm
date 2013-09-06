package Genome::Model::Tools::Synthesizer::SizeBam;

use strict;
use warnings;

use Genome;

class Genome::Model::Tools::Synthesizer::SizeBam {
    is => 'Genome::Model::Tools::Synthesizer::Base',
    has_input => [
        bam_file => {
            is  => 'String',
            doc => 'Input BAM file of alignments.',
        },
        filtered_bam_file => {
            is  => 'String',
            doc => 'Output BAM file of filtered read alignments .',
            is_output =>1
        },

        read_size_bin => {
            is => 'Text',
            doc =>'Min_max read length bin separated by \'_\' : eg 17_75',
        },

    ],
};

sub help_brief {
    "Run the Synthesizer SizeBam module to fractionate and filter bam based on specifed read lengths ";
}

sub help_detail {
    "SizeBam module of Synthesizer filters and splits the main alignment BAM into smaller size fractions based on read-lengths.Use this module to filter based on size of species. For longer reads upto 75bp, we recommend using the following bins: 17_25,26_35,36_60,61_70,71_75,17_75. This module also calculates samtools flagstat metrics for each of the filtered & sized bams";
}


sub execute {
    my $self = shift;

    my $bam_file     = $self->bam_file;
    my $filtered_bam = $self->filtered_bam_file;

    ############## FILTERING AND SIZE-SELECTING BAM ####################

    if ( defined( $self->read_size_bin ) ) {
        my $bin = $self->read_size_bin;
        my $sized_bam = filter_bam ( $bam_file, $filtered_bam,$bin);

        ############## INDEXING AND CALCULATING FLAGSTAT ###############
        unless (-s $sized_bam) {
            die( 'Failed to open filtered BAM file: ' . $sized_bam );
        }

        my $index_file= $sized_bam.'.bai';
        my $index_cmd = 'samtools index ' . $sized_bam .' > '. $index_file;
        Genome::Sys->shellcmd(
            cmd => $index_cmd,
            input_files => [$sized_bam],
            output_files => [$index_file],
        );

        my $flagstat_file= $sized_bam.'.flagstat';
        my $cmd = 'samtools flagstat '. $sized_bam .' > '. $flagstat_file ;
        Genome::Sys->shellcmd(
            cmd => $cmd,
            input_files => [$sized_bam],
            output_files => [$flagstat_file],
        );
    }
    else {
        print  "Please define a min_max size bin to filter BAM ; it should be <minsize>_<maxsize>";
    }

    return 1;
}

sub filter_bam {

    my ($bam_file, $new_bam,$bin) = @_;

    my $in_bam = Bio::DB::Bam->open($bam_file);
    unless ($in_bam) {
        die( 'Failed to open BAM file: ' . $bam_file );
    }

    my $out_bam = Bio::DB::Bam->open( $new_bam, 'w' );
    unless ($out_bam) {
        die( 'Failed to open output BAM file: ' . $new_bam );
    }

    my $header = $in_bam->header;
    $out_bam->header_write($header);

    if ( $bin =~ m/^(\d+)(_)(\d+)/ ) {
        my ( $min, $max ) = split( '_', $bin );
        while ( my $in_read = $in_bam->read1 ) {
            my $read_length = $in_read->l_qseq;
            my $in_flag     = $in_read->flag;
            my $map_score   = $in_read->qual;

            if (   $read_length > ( $min - 1 ) && $read_length < ( $max + 1 ) && $in_flag != 4 ){
                if ( $map_score eq '0' ) {
                    if ( $in_read->aux_get("XA") ) {
                        $out_bam->write1($in_read);
                    }
                }
                else {
                    $out_bam->write1($in_read);}
            }
        }
        return $new_bam;
    }
    else{
        print  "Please check format of the min_max size bin ; it should be <minsize>_<maxsize>";
    }
};

1;
