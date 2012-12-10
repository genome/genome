use strict;
use warnings;

use above 'Genome';
use Genome::Model::Tools::Music::T qw(run_test_case);

my $input_dir = Genome::Model::Tools::Music::T->input_dir;
my $actual_output_dir = Genome::Model::Tools::Music::T->output_dir;
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
