#!/usr/bin/env genome-perl

use above 'Genome';
use Test::More;
use File::Basename qw/dirname/;

if (Genome::Config->arch_os ne 'x86_64') {
    plan skip_all => 'requires 64-bit machine';
}
else {
    plan tests => 4;
}

my $pkg = 'Genome::Model::Tools::Joinx::VcfMergeForBackfill';
use_ok($pkg);

my $tmpdir = File::Temp::tempdir('joinx-VcfMergeForBackfill-XXXXX', DIR => "$ENV{GENOME_TEST_TEMP}/", CLEANUP => 1);
my $output_file = $tmpdir."/snvs.merged.vcf.gz";
my $base_dir = $ENV{GENOME_TEST_INPUTS} . "/Genome-Model-Tools-Joinx-VcfMergeForBackfill";
my $expected = $base_dir."/expected/snvs.merged.vcf.gz";
my @input_files = glob($base_dir."/inputs/*");

my $cmd = $pkg->create(
    input_files => [ @input_files ],
    output_file => $output_file,
    use_bgzip => 1,
);

ok($cmd, 'Created command');
ok($cmd->execute, 'Executed command');

print "$expected\n";
my $diff = Genome::Sys->diff_file_vs_file($output_file, $expected);
ok(!$diff, 'output matched expected result') or diag("diff results:\n" . $diff);

done_testing();
