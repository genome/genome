package Genome::Model::SomaticCapture::Report::FileSummary;

use strict;
use warnings;

use Genome;

class Genome::Model::SomaticCapture::Report::FileSummary {
    is => 'Genome::Model::Somatic::Report::FileSummary', #If more pipelines start using this, should extract to an abstract parent report class.
};

sub files_to_report {
qw(
sniper_snp_output_adaptor
sniper_snp_output_filter
sniper_snp_output_filter_hc
sniper_snp_output_filter_hc_somatic
sniper_snp_output_filter_hc_loh

adaptor_output_indel
filter_indel_output

varscan_adaptor_snp
varscan_snp_germline
varscan_snp_loh
varscan_snp_somatic

varscan_adaptor_indel
varscan_indel_germline
varscan_indel_loh
varscan_indel_somatic

merged_snp_output
merged_indel_output

merged_snp_output_filter
merged_snp_output_filter_fail
merged_snp_output_novel

annotate_output_snp
ucsc_output
ucsc_unannotated_output

annotate_output_indel

tier_1_snp_file
tier_2_snp_file
tier_3_snp_file
tier_4_snp_file
tier_1_indel_file

tier_1_snp_file_high
tier_1_snp_file_highest
tier_1_indel_file_high
tier_1_indel_file_highest
);
}

1;

