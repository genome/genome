use strict;
use warnings;

use above "Genome";
use Fcntl ':mode';
use Genome::Utility::Test qw(abort run_ok);
use List::Util qw(first);
use POSIX qw(getgroups);

use Test::More tests => 10;

eval {
    my $umask = umask;

    ok(getgrnam($ENV{GENOME_SYS_GROUP}), 'GENOME_SYS_GROUP is set to existing group') or abort;

    my $sys_group = first { $_ eq $ENV{GENOME_SYS_GROUP} } get_group_names();
    ok($sys_group, "user belongs to GENOME_SYS_GROUP ($sys_group)") or abort;

    my $test_group = first { $_ ne $sys_group } get_group_names();
    ok($test_group, "user belogs to some other group besides GENOME_SYS_GROUP ($sys_group)");

    # setup
    my $td_path = File::Temp->newdir();
    ok(-d $td_path, 'made a temp directory to work in') or abort;
    run_ok(['chown', ":$test_group", $td_path], "chown that directory to :$test_group") or abort;
    run_ok(['chmod', 2775, $td_path], 'chmod that directory to have guid sticky') or abort;

    # make sure create_directory overrides umask
    my $cd_path = File::Spec->join($td_path, 'cd');
    Genome::Sys->create_directory($cd_path);
    ok(-d $cd_path, 'made a subdirectory');
    is(gid($cd_path), Genome::Sys::gidgrnam($sys_group), 'created directory has expected gid');

    # verify mkdir, without overrides create_directory has, does not
    my $mkdir_path = File::Spec->join($td_path, 'mkdir');
    mkdir $mkdir_path;
    ok(-d $mkdir_path, 'made a subdirectory');
    is(gid($mkdir_path), Genome::Sys::gidgrnam($test_group), 'created directory has expected gid');
};

sub group_write {
    my $path = shift;
    my $mode = (stat($path))[2];
    my $perms = S_IMODE($mode);
    my $group_write = ($perms & S_IWGRP) >> 3;
    return $group_write;
}

sub gid {
    my $path = shift;
    return (stat($path))[5];
}

sub get_group_names {
    return map { (getgrgid($_))[0] } getgroups();
}
