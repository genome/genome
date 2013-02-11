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

my $pkg = 'Genome::Model::Tools::Joinx::Sort';
use_ok($pkg);

my $input = __FILE__ . ".input";
my $expected = __FILE__ . ".expected";
my $tmpdir = File::Temp::tempdir('joinx-sort-XXXXX', CLEANUP => 1, TMPDIR => 1);
my $output = join('/', $tmpdir, 'output');

my $cmd = $pkg->create(
    input_files => [$input],
    output_file => $output
    );

ok($cmd, 'Created command');
ok($cmd->execute, 'Executed command');

my $diff = Genome::Sys->diff_file_vs_file($output, $expected);
ok(!$diff, 'output matched expected result') or diag("diff results:\n" . $diff);

done_testing();
