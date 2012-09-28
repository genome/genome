use strict;
use warnings;

use above 'Genome';
use Genome::Model::Tools::Music::T qw(run_test_case $input_dir $actual_output_dir);

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
