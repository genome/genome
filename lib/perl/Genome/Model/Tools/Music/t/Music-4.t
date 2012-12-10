use strict;
use warnings;

use above 'Genome';
use Genome::Model::Tools::Music::T qw(run_test_case);

my $input_dir = Genome::Model::Tools::Music::T->input_dir;
my $actual_output_dir = Genome::Model::Tools::Music::T->output_dir;
run_test_case(
    run => "music proximity\n"
         . " --maf-file $input_dir/ucec_test.maf\n"
         . " --output-dir $actual_output_dir\n"
         . " --noskip-non-coding",
    expect => [
        'proximity_report'
    ],
);
