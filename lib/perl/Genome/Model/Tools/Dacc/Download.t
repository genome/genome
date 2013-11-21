#!/usr/bin/env genome-perl

use strict;
use warnings;

use above 'Genome';

use Test::More;

use_ok('Genome::Model::Tools::Dacc::Download') or die;

# setup
my $dir = $ENV{GENOME_TEST_INPUTS} . '/Genome-Model-Tools-Dacc';
my (%files, $expected_cmd);
my $is_running_in_lsf = 0;
no warnings qw/ once redefine /;
*Genome::Model::Tools::Dacc::available_files_and_sizes = sub{ return %files; };
*Genome::Model::Tools::Dacc::is_running_in_lsf = sub{ return $is_running_in_lsf; };
*Genome::Model::Tools::Dacc::is_host_a_blade = sub{ return 1; };
*Genome::Sys::shellcmd = sub{
    my ($class, %params) = @_;
    diag($params{cmd});
    is(
        $params{cmd}, 
        $expected_cmd,
        'Expected command matches',
    );
    return 1;
};
use warnings;

# Fail - no files on DACC
diag('Fail b/c no files on DACC');
my $dl = Genome::Model::Tools::Dacc::Download->create(
    dacc_directory => 'DACC_DIR',
    destination => $dir,
    files => [qw/ a b /],
);
ok($dl, 'create');
$dl->dump_status_messages(1);
ok(!$dl->execute, 'execute: failed as expected');

diag('Fail b/c file "b" does not exist on the DACC');
$files{a} = 2;
$dl = Genome::Model::Tools::Dacc::Download->create(
    dacc_directory => 'DACC_DIR',
    destination => $dir,
    files => [qw/ a b /],
);
ok($dl, 'create');
$dl->dump_status_messages(1);
ok(!$dl->execute, 'execute: failed as expected');

diag('Success: launch to LSF');
$files{b} = 2;
$expected_cmd = "bsub -q $ENV{GENOME_LSF_QUEUE_BUILD_WORKER} -u " . Genome::Config->user_email . " -R 'rusage[internet_download_mbps=100]' gmt dacc download /DACC_DIR/ $dir a b";
$dl = Genome::Model::Tools::Dacc::Download->create(
    dacc_directory => 'DACC_DIR',
    destination => $dir,
    files => [qw/ a b /],
    launch_to_lsf => 1,
);
ok($dl, 'create');
$dl->dump_status_messages(1);
ok($dl->execute, 'execute');

diag('Success: download');
$is_running_in_lsf = 1;
$expected_cmd = 'ascp -Q -l100M -i /gsc/scripts/share/certs/dacc/dacc.ppk jmartin@aspera.hmpdacc.org:/DACC_DIR/a jmartin@aspera.hmpdacc.org:/DACC_DIR/b '.$dir;
$dl = Genome::Model::Tools::Dacc::Download->create(
    dacc_directory => 'DACC_DIR',
    destination => $dir,
    files => [qw/ a b /],
);
ok($dl, 'create');
$dl->dump_status_messages(1);
ok($dl->execute, 'execute');

diag('Fail: download b/c file size is different');
$files{b} = 1;
$dl = Genome::Model::Tools::Dacc::Download->create(
    dacc_directory => 'DACC_DIR',
    destination => $dir,
    files => [qw/ a b /],
);
ok($dl, 'create');
$dl->dump_status_messages(1);
ok($dl->execute, 'execute');

done_testing();


=pod

=head1 Tests

=head1 Disclaimer

 Copyright (C) 2010 Washington University Genome Sequencing Center

 This script is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY or the implied warranty of MERCHANTABILITY
 or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public
 License for more details.

=head1 Author(s)

 Eddie Belter <ebelter@watson.wustl.edu>

=cut
