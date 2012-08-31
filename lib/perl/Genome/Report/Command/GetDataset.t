#!/usr/bin/env genome-perl

use strict;
use warnings;

use above 'Genome';

use Test::More;

use_ok('Genome::Report::Command::GetDataset') or die;

my @output_types = Genome::Report::Command::GetDataset->output_types;
for my $output_type ( @output_types ) {
    my $cmd = Genome::Report::Command::GetDataset->create(
        report_directory => $ENV{GENOME_TEST_INPUTS} . '/Genome-Report-XSLT/Assembly_Stats',
        dataset_name => 'stats',
        output_type => $output_type,
    );
    ok($cmd, 'create get dataset command for '.$output_type);
    $cmd->dump_status_messages(1);
    ok($cmd->execute, 'execute');
}

done_testing();
exit;

