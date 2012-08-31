#!/usr/bin/env genome-perl

use strict;
use warnings;

use above 'Genome';

use Test::More;

use_ok('Genome::Report::Command::Xslt');

no warnings;
local *Genome::Report::XSLT::transform_report = sub { # so we don't test this twice
    my $content = <<EOS;

Summary for Amplicon Assembly Model (Name:  Build Id:)

------------------------
Stats
------------------------

Attempted
Assembled 5
Assembly Success 100.00%

Length Average 1399
Length Median 1396
Length Maximum 1413
Length Minimum 1385

Quality Base Average 62.75
Quality >= 20 per Assembly 1349.80

Reads Assembled 20
Reads Total 30
Reads Assembled Success 66.67%
Reads Assembled Average 4.00
Reads Assembled Median 3
Reads Assembled Maximum 6
Reads Assembled Minimum 3

------------------------

For full report, including quality hisotgram go to:
http://

EOS
    return { 
        media_type => 'text/plain',
        output_type => 'txt',
        content => $content,
    };
};
use warnings;

my $dir = $ENV{GENOME_TEST_INPUTS} . '/Genome-Report-XSLT';
my $tmpdir = Genome::Sys->base_temp_directory;
my $cmd = Genome::Report::Command::Xslt->create(
    report_directory => $dir.'/Assembly_Stats',
    xsl_file => $dir.'/AssemblyStats.txt.xsl',
    output_file => $tmpdir.'/Assembly_Stats.txt',
);
ok($cmd, 'create list dataset command');
$cmd->dump_status_messages(1);
ok($cmd->execute, 'execute');
ok(-e $cmd->output_file, 'Output file exists');

done_testing();
exit;

