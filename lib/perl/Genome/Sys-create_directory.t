use strict;
use warnings;

use above "Genome";
use Fcntl ':mode';
use Genome::Utility::Test qw(abort run_ok);

use Test::More tests => 8;

eval {
    my $umask = umask;

    my $test_group = 'apipe';

    # setup
    my $td_path = File::Temp->newdir();
    ok(-d $td_path, 'made a temp directory to work in') or abort;
    isnt($ENV{GENOME_SYS_GROUP}, $test_group, "test assumes GENOME_SYS_GROUP is not set to $test_group") or abort;
    run_ok(['chown', ":$test_group", $td_path], "chown that directory to :$test_group") or abort;
    run_ok(['chmod', 2775, $td_path], 'chmod that directory to have guid sticky') or abort;

    # make sure create_directory overrides umask
    my $cd_path = File::Spec->join($td_path, 'cd');
    Genome::Sys->create_directory($cd_path);
    ok(-d $cd_path, 'made a subdirectory');
    is(gid($cd_path), Genome::Sys::gidgrnam($ENV{GENOME_SYS_GROUP}), 'created directory has expected gid');

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
