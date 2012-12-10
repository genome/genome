use strict;
use warnings;

use above 'Genome';
use Genome::Model::Tools::Music::T qw(run_test_case);

my $input_dir = Genome::Model::Tools::Music::T->input_dir;
my $actual_output_dir = Genome::Model::Tools::Music::T->output_dir;
run_test_case(
    skip => q(Skip this tool, because it's crap and needs to be rewritten.),
    run => "music cosmic-omim \n"
         . " --maf-file $input_dir/ucec_test.maf\n"
         . " --output-file $actual_output_dir/cosmic_omim.maf\n"
         . " --no-verbose",
    expect => [
        'cosmic_omim.maf'
    ],
);
