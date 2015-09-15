package Genome::Qc::Config;

use strict;
use warnings;
use Genome;

class Genome::Qc::Config {
    is => 'UR::Value',
    id_by => [
        name => {
            is => 'String',
        },
    ],
};

sub get_commands_for_alignment_result {
    my $self = shift;
    my $is_capture = shift;

    my %config = (
        picard_collect_gc_bias_metrics => {
            class => 'Genome::Qc::Tool::Picard::CollectGcBiasMetrics',
            params => {
                input_file => 'bam_file',
                refseq_file => 'reference_sequence',
                assume_sorted => 1,
                use_version => 1.123,
                output_file=> 'output_file',
                chart_output => 'chart_output',
                # metric_accumulation_level => ['SAMPLE'], #not supported in 1.123
            },
        },
        picard_mark_duplicates => {
            class => 'Genome::Qc::Tool::Picard::MarkDuplicates',
            params => {
                output_file => 'output_file',
                input_file => 'bam_file',
                use_version => 1.123,
            },
        },
        samtools_flagstat => {
            class => 'Genome::Qc::Tool::Samtools::Flagstat',
            params => {
                'bam-file' => 'bam_file',
            },
        },
        picard_collect_alignment_summary_metrics => {
            class => 'Genome::Qc::Tool::Picard::CollectAlignmentSummaryMetrics',
            params => {
                input_file => 'bam_file',
                refseq_file => 'reference_sequence',
                use_version => 1.123,
                metric_accumulation_level => ['SAMPLE', 'READ_GROUP'],
            }
        },
        picard_collect_insert_size_metrics => {
            class => 'Genome::Qc::Tool::Picard::CollectInsertSizeMetrics',
            params => {
                input_file => 'bam_file',
                histogram_file => 'histogram_file',
                use_version => 1.123,
                metric_accumulation_level => ['SAMPLE'],
            }
        },
    );

    if ($is_capture) {
        $config{picard_calculate_hs_metrics} = {
            class => 'Genome::Qc::Tool::Picard::CalculateHsMetrics',
            params => {
                input_file => 'bam_file',
                bait_intervals => 'bait_intervals', #region_of_interest_set
                target_intervals => 'target_intervals', #target_region_set
                use_version => 1.123,
                metric_accumulation_level => ['SAMPLE'],
            },
        };
    }
    else {
        $config{picard_collect_wgs_metrics} = {
            class => 'Genome::Qc::Tool::Picard::CollectWgsMetrics',
            params => {
                input_file => 'bam_file',
                reference_sequence => 'reference_sequence',
                use_version => 1.123,
            },
        };
        #Ultimately we also want to run verifyBamId on exome data but we need
        #to preprocess the vcf to limit to ROI
        $config{verify_bam_id} = {
            class => 'Genome::Qc::Tool::VerifyBamId',
            params => {
                vcf_file => '/gscmnt/gc2802/halllab/abelhj/gatk_utah_041815/indiv/vbid/Omni25_genotypes_1525_samples_v2.b37.PASS.ALL.sites.vcf.gz',
                bam => 'bam_file',
                max_depth => '150',
                precise => '1',
                version => '20120620',
                ignore_read_group => 0,
            }
        };
    }

    return \%config;
}

1;

