#!/usr/bin/env genome-perl

use strict;
use warnings;

use above 'Genome';

use Test::More;

use_ok('Genome::Report::Command::ListDatasets') or die;

my $cmd = Genome::Report::Command::ListDatasets->create(
    report_directory => $ENV{GENOME_TEST_INPUTS} . '/Genome-Report-XSLT/Assembly_Stats',
);
ok($cmd, 'create list dataset command');
$cmd->dump_status_messages(1);
ok($cmd->execute, 'execute');

done_testing();
