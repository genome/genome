#!/usr/bin/env genome-perl

use strict;
use warnings;

use above 'Genome';

use Test::More;

use_ok('Genome::Model::Tools::Dacc::Upload') or die;

# setup
my $dir = $ENV{GENOME_TEST_INPUTS} . '/Genome-Model-Tools-Dacc';
my @files_to_upload = map { $dir.'/'.$_ } (qw/ a b /);
my %files;
my $expected_cmd;
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

diag('Fail b/c files do not exist');
my $up = Genome::Model::Tools::Dacc::Upload->create(
    dacc_directory => 'DACC_DIR',
    files => [qw/ DOES_NOT_EXIST /],
);
ok($up, 'create');
$up->dump_status_messages(1);
ok(!$up->execute, 'execute: failed as expected');

diag('Success: launch to LSF');
$expected_cmd = "bsub -q $ENV{GENOME_LSF_QUEUE_BUILD_WORKER} -u " . Genome::Config->user_email . " -R 'rusage[internet_upload_mbps=100,aspera_upload_mbps=100]' gmt dacc upload /DACC_DIR/ $dir/a $dir/b";
$up = Genome::Model::Tools::Dacc::Upload->create(
    dacc_directory => 'DACC_DIR',
    files => \@files_to_upload,
    launch_to_lsf => 1,
);
ok($up, 'create');
$up->dump_status_messages(1);
ok($up->execute, 'execute');

diag('Success: upload');
%files = ( a => 2, b => 2 );
$is_running_in_lsf = 1;
$expected_cmd = 'ascp -Q -l100M -i /gsc/scripts/share/certs/dacc/dacc.ppk -d '.join(' ', @files_to_upload).' jmartin@aspera.hmpdacc.org:/DACC_DIR/';
$up = Genome::Model::Tools::Dacc::Upload->create(
    dacc_directory => 'DACC_DIR',
    files => \@files_to_upload,
);
ok($up, 'create');
$up->dump_status_messages(1);
ok($up->execute, 'execute');

diag('Fail: upload b/c file size is different');
$files{b} = 1;
$up = Genome::Model::Tools::Dacc::Upload->create(
    dacc_directory => 'DACC_DIR',
    files => \@files_to_upload,
);
ok($up, 'create');
$up->dump_status_messages(1);
ok($up->execute, 'execute');

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
