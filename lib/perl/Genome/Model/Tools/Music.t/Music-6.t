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
    run => "music smg\n"
         . " --gene-mr-file $input_dir/gene_mrs\n"
         . " --output-file $actual_output_dir/smgs\n"
         . " --max-fdr 0.55",
    expect => [
        'smgs',
        'smgs_detailed'
    ],
);
