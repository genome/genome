use strict;
use warnings;

use above 'Genome';
use Cwd qw(realpath);
use File::Basename qw(dirname);
use lib realpath(dirname(__FILE__));
use Test::Music qw(run_test_case);

my $input_dir = Test::Music->input_dir;
my $actual_output_dir = Test::Music->output_dir;
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
