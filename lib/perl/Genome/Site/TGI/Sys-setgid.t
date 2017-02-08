use strict;
use warnings;

use above "Genome";
use Genome::Utility::File::Mode qw(mode);

use autodie qw(chown);
use Fcntl ':mode';
use File::Spec qw();
use List::Util qw(first);
use POSIX qw(getgroups);
use Sub::Install qw(reinstall_sub);

use File::Spec;

use Test::More tests => 3;

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

my $nfs_working_dir = File::Spec->join(
    Genome::Config::get('site_lock_dir'),
    'temp_for_setgid_test'
);


subtest 'create_directory (on NFS) preserves setgid' => sub {
    plan tests => 2;

    my $dir = File::Spec->join(
        Genome::Config::get('site_lock_dir'),
        'temp_for_setgid_test'
    );

    my $base_dir = File::Temp->newdir(DIR => $nfs_working_dir);
    my $setgid_dirname = File::Spec->join($base_dir->dirname, 'setgid');
    Genome::Sys->create_directory($setgid_dirname);
    my $mode = mode($setgid_dirname);
    $mode->add_setgid();
    ok(-g $setgid_dirname, 'setgid_dirname has setgid set');

    my $new_dirname = File::Spec->join($setgid_dirname, 'genome');
    Genome::Sys->create_directory($new_dirname);
    ok(-g $new_dirname, 'new_dirname has setgid set');
};

subtest 'create_directory (on NFS) preserves setgid even when set_gid is needed' => sub {
    plan tests => 6;

    my $set_gid = \&Genome::Sys::set_gid;
    my $set_gid_ran = 0;
    local *Genome::Sys::set_gid;
    reinstall_sub({
        into => 'Genome::Sys',
        as => 'set_gid',
        code => sub { $set_gid->(@_); $set_gid_ran++ },
    });

    ok(getgrnam(Genome::Config::get('sys_group')), 'sys_group is set to existing group') or return;

    my $sys_group = first { $_ eq Genome::Config::get('sys_group') } get_group_names();
    ok($sys_group, "user belongs to sys_group ($sys_group)") or return;

    my $test_group = first { $_ ne $sys_group } get_group_names();
    ok($test_group, "user belongs to some other group besides sys_group ($sys_group)") or return;

    my $base_dir = File::Temp->newdir(DIR => $nfs_working_dir);
    my $setgid_dirname = File::Spec->join($base_dir->dirname, 'setgid');
    Genome::Sys->create_directory($setgid_dirname);
    chown -1, gidgrnam($test_group), $setgid_dirname;
    my $mode = mode($setgid_dirname);
    $mode->add_setgid();
    ok(-g $setgid_dirname, 'setgid_dirname has setgid set');

    my $new_dirname = File::Spec->join($setgid_dirname, 'genome');
    $set_gid_ran = 0;
    Genome::Sys->create_directory($new_dirname);
    is($set_gid_ran, 1, 'set_gid ran');
    ok(-g $new_dirname, 'new_dirname has setgid set');
};

sub get_group_names {
    return map { (getgrgid($_))[0] } getgroups();
}

sub gidgrnam {
    my $group_name = shift;
    (getgrnam($group_name))[2];
}
