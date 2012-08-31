#!/usr/bin/env genome-perl

use strict;
use warnings;

use above 'Genome';

use Test::More;

no warnings;
# overload 'Close' to not send the mail, but to cancel it 
use Mail::Sender;
*Mail::Sender::Close = sub{ my $sender = shift; $sender->Cancel; return 1; };
use warnings;

use_ok('Genome::Report::Command::Email') or die;

my $dir = $ENV{GENOME_TEST_INPUTS} . '/Genome-Report-XSLT';
my $cmd = Genome::Report::Command::Email->create(
    report_directory => $dir.'/Assembly_Stats',
    xsl_files => $dir.'/AssemblyStats.txt.xsl',
    to => Genome::Config->user_email,
);
ok($cmd, 'create email command');
$cmd->dump_status_messages(1);
ok($cmd->execute, 'execute');

done_testing();
exit;

