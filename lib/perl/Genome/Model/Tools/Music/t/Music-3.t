use strict;
use warnings;

use above 'Genome';
use Genome::Model::Tools::Music::T qw(run_test_case $input_dir $actual_output_dir);

run_test_case(
    run => "music pfam\n"
         . " --maf-file $input_dir/ucec_test.maf\n"
         . " --output-file $actual_output_dir/pfam.maf",
    expect => [
        'pfam.maf'
    ],
);
