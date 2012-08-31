#!/usr/bin/env genome-perl

use above 'Genome';
use Test::More;

if (Genome::Config->arch_os ne 'x86_64') {
    plan skip_all => 'requires 64-bit machine';
}
else {
    plan tests => 4;
}

my $pkg = 'Genome::Model::GenotypeMicroarray::Command::CreateGoldSnpBed';
use_ok($pkg);

my $input = __FILE__ . ".input";
my $expected = __FILE__ . ".expected";
my $tmpdir = File::Temp::tempdir('create-gold-snp-bed-XXXXX', DIR => "$ENV{GENOME_TEST_TEMP}/", CLEANUP => 1);
my $output = join('/', $tmpdir, 'output');

my $ref = Genome::Model::Build::ImportedReferenceSequence->get(
    name => "NCBI-human-build36"
    );

my $cmd = $pkg->create(
    input_file => $input,
    output_file => $output,
    reference => $ref,
    );

ok($cmd, 'Created command');
ok($cmd->execute, 'Executed command');

my $diff = Genome::Sys->diff_file_vs_file($output, $expected);
ok(!$diff, 'output matched expected result') or diag("diff results:\n" . $diff);

done_testing();
