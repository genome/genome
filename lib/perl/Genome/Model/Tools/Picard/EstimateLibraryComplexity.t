
#!/usr/bin/env genome-perl

use strict;
use warnings;

use above 'Genome';
use Genome::Utility::Test 'compare_ok';
use Test::More;
use File::Spec qw();

my $pkg = 'Genome::Model::Tools::Picard::EstimateLibraryComplexity';
use_ok($pkg);

my $data_dir = sprintf "%s.d", __FILE__;

my $input_file = File::Spec->catfile($data_dir, "simulated.sam");
my $expected_file = File::Spec->catfile($data_dir, "expected.txt");
my $output_file = Genome::Sys->create_temp_file_path;

my $cmd = $pkg->create(
    input_file => $input_file,
    output_file => $output_file,
    min_mean_quality => 1,
    );

ok($cmd, "created command");
ok($cmd->execute, "executed command");

compare_ok($expected_file, $output_file, filters => ['^#.*']);

done_testing();
