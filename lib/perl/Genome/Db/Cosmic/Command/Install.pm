#!/usr/bin/env genome-perl
use strict;
use warnings;

package Genome::Db::Cosmic::Command::Install;

class Genome::Db::Cosmic::Command::Install {
    is => 'Command::V2',
    has => [
        version => { is => 'Text',
                    shell_args_position => 1,
                    doc => 'the version of cosmic to install, corresponding to a branch in the genome-db-cosmic-data repository' },
        repo    => { is => 'Text',
                    is_optional => 1,
                    default_value => 'git@github.com:genome-vendor/genome-db-cosmic-data.git',
                    doc => 'override the repository URL' },
    ],
    doc => 'install some version of COSMIC that has been imported into the GMS, and for which there is a branch in the canonical repository'
};

sub execute {
    my $self = shift;
    my @dirs = split(":",$ENV{GENOME_DB});
    unless (@dirs) {
        die "The GENOME_ENV environment variable must be set to the location of file-based databases!";
    }
    $self->status_message("Found database directories: @dirs");

    my $dir = $dirs[0];
    $self->status_message("Selected database directory: $dir");

    my $repo = $self->repo;
    my $version = $self->version;

    my $subdir = $version;
    $subdir =~ s/\./_v/;

    my $cmd = "cd $dir; git clone $repo -b $version $subdir";
    Genome::Sys->shellcmd(cmd => $cmd);

    my $latest = $dir . '/latest';
    if (-e $latest) {
        unlink $latest;
    }
    symlink "$dir/$subdir", $latest;

    return 1;
}

1;
