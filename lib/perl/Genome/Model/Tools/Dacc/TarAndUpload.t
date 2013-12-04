#!/usr/bin/env genome-perl

use strict;
use warnings;

use above 'Genome';

use Test::More;

use_ok('Genome::Model::Tools::Dacc::TarAndUpload') or die;

my $dir = $ENV{GENOME_TEST_INPUTS} . '/Genome-Model-Tools-Dacc';
my @files = map { $dir.'/'.$_ } (qw/ a b /);
my $tar_file = $dir.'/dacc.tar.gz';

my $cnt = 0;
my @cmds = (
    "tar cvzf $tar_file ".join(' ', @files),
    "bsub -q $ENV{GENOME_LSF_QUEUE_BUILD_WORKER} -u " . Genome::Config->user_email . " -R 'rusage[internet_upload_mbps=100,aspera_upload_mbps=100]' gmt dacc upload /DACC_DIR/ $tar_file"
);
no warnings qw/ once redefine /;
*Genome::Sys::shellcmd = sub{
    my ($class, %params) = @_;
    diag($params{cmd});
    is($params{cmd}, $cmds[$cnt], 'Command matches');
    $cnt++;
    return 1;
};
use warnings;

my $tnl = Genome::Model::Tools::Dacc::TarAndUpload->create(
    dacc_directory => 'DACC_DIR',
    files => \@files,
    tar_file => $tar_file,
);
ok($tnl, 'create');
$tnl->dump_status_messages(1);
ok($tnl->execute, 'execute');

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
