package Genome::Model::Tools::Cufflinks::FpkmTrackingStats;

use strict;
use warnings;

use Genome;
use Statistics::Descriptive;



my @HEADERS = qw/total OK LOWDATA FAIL touched min_coverage max_coverage mean_coverage stdev_coverage min_fpkm max_fpkm mean_fpkm stdev_fpkm/;

class Genome::Model::Tools::Cufflinks::FpkmTrackingStats {
    is => ['Command'],
    has => [
        fpkm_tracking_file => {
            is => 'Text',
            doc => 'An isoforms.fpkm_tracking file output from cufflinks.',
        },
        output_summary_tsv => {
            is => 'Text',
            doc => 'An output summary file with tab-separated values.',
            is_optional => 1,
        },
        _stats_hash_ref => {
            is_optional => 1,
        }
    ],
};

sub headers {
    return @HEADERS;
}

sub execute {
    my $self = shift;

    my $reader = Genome::Utility::IO::SeparatedValueReader->create(
        input => $self->fpkm_tracking_file,
        separator => "\t",
    );
    my $coverage_stats = Statistics::Descriptive::Sparse->new();
    my $fpkm_stats = Statistics::Descriptive::Sparse->new();
    my %stats;
    while (my $data = $reader->next) {
        # Example data structure for isoforms.fpkm_tracking
        # $data = {
        #  'coverage' => '0',
        #  'tss_id' => '-',
        #  'gene_id' => 'Xkr4',
        #  'status' => 'OK',
        #  'FPKM_conf_lo' => '0',
        #  'locus' => '1:3204562-3661579',
        #  'gene_short_name' => 'Xkr4',
        #  'tracking_id' => 'ENSMUST00000070533',
        #  'class_code' => '-',
        #  'length' => '3634',
        #  'FPKM' => '0',
        #  'nearest_ref_id' => '-',
        #  'FPKM_conf_hi' => '0'
        #};
        $stats{total}++;
        my $status = $data->{status};
        $stats{$status}++;
        if ($status eq 'OK') {
            my $coverage = $data->{coverage};
            my $length = $data->{length};
            # I don't think this will actually work given the normalization and bias correction performed during cufflinks
            # my $reads = int($coverage * $length);
            if ($coverage =~ /\d+/) {
                if ($coverage > 0) {
                    $stats{touched}++;
                    # The coverage metrics are now weighted by length
                    my @data;
                    for (1 .. $length) {
                        push @data, $coverage;
                    }
                    $coverage_stats->add_data(@data);
                }
            }
            my $fpkm = $data->{FPKM};
            if ($fpkm > 0) {
                $fpkm_stats->add_data($fpkm);
            }
        }
    }
    $stats{min_coverage} = $coverage_stats->min;
    $stats{max_coverage} = $coverage_stats->max;
    $stats{mean_coverage} = $coverage_stats->mean;
    $stats{stdev_coverage} = $coverage_stats->standard_deviation;

    $stats{min_fpkm} = $fpkm_stats->min;
    $stats{max_fpkm} = $fpkm_stats->max;
    $stats{mean_fpkm} = $fpkm_stats->mean;
    $stats{stdev_fpkm} = $fpkm_stats->standard_deviation;

    my @headers = $self->headers;
    my $writer = Genome::Utility::IO::SeparatedValueWriter->create(
        output => $self->output_summary_tsv,
        separator => "\t",
        headers => \@headers,
    );
    $writer->write_one(\%stats);
    $self->_stats_hash_ref(\%stats);
    return 1;
};


1;
