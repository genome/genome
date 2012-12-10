use strict;
use warnings;

use above 'Genome';
use Genome::Model::Tools::Music::T qw(run_test_case);

my $input_dir = Genome::Model::Tools::Music::T->input_dir;
my $actual_output_dir = Genome::Model::Tools::Music::T->output_dir;
run_test_case(
    skip => q(Skip this tool, R arbitrarily changes precision causing expected output to differ.),
    run => "music mutation-relation\n"
         . " --bam-list $input_dir/bam_list\n"
         . " --maf-file $input_dir/ucec_test.maf\n"
         . " --gene-list $input_dir/smgs\n"
         . " --output-file $actual_output_dir/mutation_relations\n"
         . " --mutation-matrix-file $actual_output_dir/mutation_matrix\n"
         . " --permutations 500\n"
         . " --noskip-non-coding",
    expect => [
        'mutation_relations',
        'mutation_matrix'
    ],
);
