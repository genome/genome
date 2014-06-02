use strict;
use warnings;

use above "Genome";
use Genome::Utility::File::Mode qw(mode);

use Fcntl ':mode';
use File::Spec qw();

use Test::More tests => 2;

subtest 'create_directory preserves setgid' => sub {
    plan tests => 2;

    my $base_dir = File::Temp->newdir();
    my $setgid_dirname = File::Spec->join($base_dir->dirname, 'setgid');
    Genome::Sys->create_directory($setgid_dirname);
    my $mode = mode($setgid_dirname);
    $mode->add_setgid();
    ok(-g $setgid_dirname, 'setgid_dirname has setgid set');

    my $new_dirname = File::Spec->join($setgid_dirname, 'genome');
    Genome::Sys->create_directory($new_dirname);
    ok(-g $new_dirname, 'new_dirname has setgid set');
};

subtest 'create_directory (on NFS) preserves setgid' => sub {
    plan tests => 2;

    my $base_dir = File::Temp->newdir(DIR => '/gsc/var/tmp');
    my $setgid_dirname = File::Spec->join($base_dir->dirname, 'setgid');
    Genome::Sys->create_directory($setgid_dirname);
    my $mode = mode($setgid_dirname);
    $mode->add_setgid();
    ok(-g $setgid_dirname, 'setgid_dirname has setgid set');

    my $new_dirname = File::Spec->join($setgid_dirname, 'genome');
    Genome::Sys->create_directory($new_dirname);
    ok(-g $new_dirname, 'new_dirname has setgid set');
};
