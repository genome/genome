use strict;
use warnings;

use above 'Genome';
use Genome::Model::Tools::Music::T qw(run_test_case $input_dir $actual_output_dir);

run_test_case(
    run => "music proximity\n"
         . " --maf-file $input_dir/ucec_test.maf\n"
         . " --output-dir $actual_output_dir\n"
         . " --noskip-non-coding",
    expect => [
        'proximity_report'
    ],
);
