package Genome::Model::Tools::BioSamtools::CleanBam;

use strict;
use warnings;

use Genome;
use YAML;

class Genome::Model::Tools::BioSamtools::CleanBam {
    is => ['Genome::Model::Tools::BioSamtools'],
    has_input => [
        input_bam_file => {
            is => 'Text',
            doc => 'An input BAM format file of alignment data'
        },
        output_bam_file => {
            is => 'Text',
            doc => 'An output BAM format file of alignment data',
        },
        summary_output_file => {
            is => 'Text',
            doc => 'An output tsv file of summary metrics',
        }
    ],
};

sub execute {
    my $self = shift;
    my $summary_fh = Genome::Sys->open_file_for_writing($self->summary_output_file);
    my $in_bam = Bio::DB::Bam->open($self->input_bam_file);
    my $out_bam = Bio::DB::Bam->open($self->output_bam_file,'w');
    my $header = $in_bam->header();
    $out_bam->header_write($header);
    my $target_names = $header->target_name();
    my $target_lengths = $header->target_len();

    my %stats_summary;
    while (my $align = $in_bam->read1()) {
        my $tid = $align->tid;
        my $target_name = $target_names->[$tid];
        my $target_length = $target_lengths->[$tid];
        my $end = $align->calend;
        if ($end >= $target_length) {
            # Modifying alignment objects with Bio::DB::Sam is not currently supported
            #my $flag = $align->flag();
            # set unmapped flag
            #$flag += 4;
            #$align->flag($flag);
            $stats_summary{removed_reads}++;
            $stats_summary{removed_bp} += $align->l_qseq;
            $stats_summary{overhanging_bp} += ($end - $target_length);
            next;
        }
        $out_bam->write1($align);
    }
    print $summary_fh Dump %stats_summary;
    return 1;
}

1;
