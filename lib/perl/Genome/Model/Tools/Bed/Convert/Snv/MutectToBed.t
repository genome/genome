#!/usr/bin/env genome-perl

use strict;
use warnings;

use Test::More;

if (Genome::Config->arch_os ne 'x86_64') {
    plan skip_all => 'requires 64-bit machine';
}

use above 'Genome';

use_ok('Genome::Model::Tools::Bed::Convert::Snv::MutectToBed');

#my $test_dir = "/gscmnt/sata831/info/medseq/dlarson/mutect_testing/Genome-Model-Tools-Bed-Convert-Snv-MutectToBed";
my $test_dir = $ENV{GENOME_TEST_INPUTS} . "/Genome-Model-Tools-Bed-Convert-Snv-MutectToBed";

my $expected_base = "expected.v2";
my $expected_dir = "$test_dir/$expected_base";
my $all_expected_file = "$expected_dir/output.all.bed";
my $hq_expected_file = "$expected_dir/output.hq.bed";
my $lq_expected_file = "$expected_dir/output.lq.bed";

my $input_file = "$test_dir/snvs.v1.1.4.hq"; #updated to same as v1.1.4 output to check results look good.

diag("Testing without any filtering");
test_conversion($all_expected_file, source => $input_file);
diag("Testing filtering to high quality only");
test_conversion($hq_expected_file, source => $input_file, limit_variants_to => 'hq');
diag("Testing filtering to low quality only");
test_conversion($lq_expected_file, source => $input_file, limit_variants_to => 'lq');
done_testing();

sub test_conversion {
    my ($expected_file, %params) = @_;
    my $output_file = Genome::Sys->create_temp_file_path;
    $params{output} = $output_file;
    my $command = Genome::Model::Tools::Bed::Convert::Snv::MutectToBed->create(%params);

    ok($command, 'Command created');
    my $rv = $command->execute;
    ok($rv, 'Command completed successfully');
    ok(-s $output_file, "output file created");

    my $diff = Genome::Sys->diff_file_vs_file($output_file, $expected_file);
    ok(!$diff, 'output matched expected result')
        or diag("diff results:\n" . $diff);
}

