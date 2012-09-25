use strict;
use warnings;

use above 'Genome';
use Genome::Model::Tools::Music::T qw(run_test_case $input_dir $actual_output_dir);

run_test_case(
    run => "music path-scan\n"
         . " --bam-list $input_dir/bam_list\n"
         . " --maf-file $input_dir/ucec_test.maf\n"
         . " --pathway-file $input_dir/kegg_db_120910\n"
         . " --gene-covg-dir $input_dir/gene_covgs\n"
         . " --output-file $actual_output_dir/sm_pathways\n"
         . " --bmr 1.8E-06",
    expect => [
        'sm_pathways',
        'sm_pathways_detailed'
    ],
);
