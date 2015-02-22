#!/usr/bin/env genome-perl

use strict;
use warnings;

use Test::More;

use File::Spec qw();
use File::Temp qw();
use File::stat qw(stat);
use POSIX qw(getgid getgroups);

use Genome;
use Genome::Sys qw();
use Genome::Utility::File::Mode qw(mode);

my $gid = getgid();
my $alt_gid = (getgroups())[1]; # 0 = getgid
unless ($alt_gid) {
    die 'no alternate gid available';
}

# The CORE::rename case demonstrates the undesired behavior compared to
# the Genome::Sys::rename case.
my @cases = (
    ['CORE::rename'       , sub { CORE::rename($_[0], $_[1]) }, $gid    , 0],
    ['Genome::Sys::rename', sub { Genome::Sys->rename(@_)    }, $alt_gid, 1],
);
plan tests => scalar(@cases);
for my $case (@cases) {
    my ($name, $fn, $exp_gid, $setgid) = @{$case};
    subtest $name => sub {
        plan tests => 3;

        my $tmp_dir = File::Temp->newdir();
        my $setgid_dir = File::Temp->newdir();
        subtest 'setup' => sub {
            plan tests => 4;

            is(stat($tmp_dir)->gid, $gid, 'tmp_dir has expected gid');
            ok(! -g $tmp_dir, 'tmp_dir does not have setgid bit');

            chown -1, $alt_gid, $setgid_dir;
            is(stat($setgid_dir)->gid, $alt_gid, 'setgid_dir has expected gid');

            mode($setgid_dir->dirname)->add_setgid();
            ok(-g $setgid_dir, 'setgid_dir has setgid bit');
        };

        subtest 'file behavior' => sub {
            plan tests => 5;

            my $tmp_filename = File::Spec->join($tmp_dir->dirname, 'somefile');
            ok(! -e $tmp_filename, 'tmp_filename does not already exist');

            Genome::Sys->touch($tmp_filename);
            ok(-f $tmp_filename, 'tmp_filename created');

            my $setgid_filename = File::Spec->join($setgid_dir->dirname, 'somefile');
            ok(! -e $setgid_filename, 'setgid_filename does not already exist');

            $fn->($tmp_filename, $setgid_filename);
            ok(-f $setgid_filename, 'file renamed') or return;

            is(stat($setgid_filename)->gid, $exp_gid, 'renamed file has expected gid');
        };

        subtest 'directory behavior' => sub {
            plan tests => 6;

            my $tmp_subdirname = File::Spec->join($tmp_dir->dirname, 'somedir');
            ok(! -e $tmp_subdirname, 'tmp_subdirname does not already exist');

            mkdir $tmp_subdirname;
            ok(-d $tmp_subdirname, 'tmp_subdirname created');

            my $setgid_subdirname = File::Spec->join($setgid_dir->dirname, 'somedir');
            ok(! -e $setgid_subdirname, 'setgid_subdirname does not already exist');

            $fn->($tmp_subdirname, $setgid_subdirname);
            ok(-d $setgid_subdirname, 'subdir renamed') or return;

            is(stat($setgid_subdirname)->gid, $exp_gid, 'renamed subdir has expected gid');
            ok((-g $setgid_subdirname == $setgid), 'renamed subdir has expected setgid value');
        };
    };
}
