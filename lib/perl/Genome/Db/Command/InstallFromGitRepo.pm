package Genome::Db::Command::InstallFromGitRepo;
use strict;
use warnings;
use Genome;

class Genome::Db::Command::InstallFromGitRepo {
    is => 'Command::V2',
    has => [
        source  => { is => 'Text',
                    doc => 'the data source' },
        version => { is => 'Text',
                    shell_args_position => 1,
                    doc => 'the version of cosmic to install, corresponding to a branch in the genome-db-$SOURCE-data repository' },
        repo    => { is => 'Text',
                    is_optional => 1,
                    doc => 'override the repository URL' },
    ],
    doc => 'install some version data set on the local GMS for which there is a branch in the canonical repository'
};

sub help_detail {
    return <<EOS
Install a file-based database from one of the canonical set at http://github.com/genome-vendor/genome-db-\$SOURCE-data.git.
EOS
}

sub _resolve_repo {
    my $self = shift;
    my $source = $self->source;
    return "http://github.com/genome-vendor/genome-db-${source}-data.git";
}

sub execute {
    my $self = shift;
    my @dirs = split(":",$ENV{GENOME_DB});
    unless (@dirs) {
        die "The GENOME_ENV environment variable must be set to the location of file-based databases!";
    }
    $self->status_message("Found database directories: @dirs");

    my $dir = $dirs[0];
    $self->status_message("Selected database directory: $dir");

    my $repo = $self->repo || $self->_resolve_repo;
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
