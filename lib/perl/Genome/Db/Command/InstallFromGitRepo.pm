package Genome::Db::Command::InstallFromGitRepo;
use strict;
use warnings;
use Genome;

class Genome::Db::Command::InstallFromGitRepo {
    is => 'Command::V2',
    has => [
        source  => { is => 'Text',
                    doc => 'the data source' },
        subsource => { is => 'Text',
                       doc => 'the data sub-source (if applicable)' },
        species => { is => 'Text',
                     doc => 'the species of the source (if applicable, e.g. human or mouse)' },
        branch => { is => 'Text',
                    doc => 'the branch of database to install, corresponding to a branch name in the genome-db-$SOURCE-data repository' },
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
    my $subsource = $self->subsource;
    if ($subsource){
        return "http://github.com/genome-vendor/genome-db-${source}-${subsource}-data.git"
    }else{
        return "http://github.com/genome-vendor/genome-db-${source}-data.git";
    }
}

sub execute {
    my $self = shift;
    my @dirs = split(":",$ENV{GENOME_DB});
    unless (@dirs) {
        die "The GENOME_ENV environment variable must be set to the location of file-based databases!";
    }
    $self->status_message("Database will be installed at GENOME_DB: " . $ENV{GENOME_DB}); 
    $self->status_message("Found database directories: @dirs");

    my $dir = $dirs[0];
    $self->status_message("Selected database directory: $dir");

    my $repo = $self->repo || $self->_resolve_repo;
    my $source = $self->source;
    my $subsource = $self->subsource;
    my $species = $self->species;
    my $branch = $self->branch;

    $self->status_message("Repo is: " . $repo);
    $self->status_message("Source is: " . $source);
    $self->status_message("Subsource is: " . $subsource) if $self->subsource;
    $self->status_message("Species is: " . $species) if $self->species;
    $self->status_message("Branch is: " . $branch);

    my $subdir = $branch;
    $subdir = $species . "/" . $subdir if $species;
    $subdir = $subsource . "/" . $subdir if $subsource;
    $subdir = $source . "/" . $subdir;

    my $clone_cmd = "cd $dir; git clone $repo -b $branch $subdir";
    Genome::Sys->shellcmd(cmd => $clone_cmd);

    my $exit_code = `bash -c 'cd $subdir; git ls-remote --exit-code . origin/$branch &> /dev/null; echo \$\?'`;
    chomp($exit_code);
    
    unless ($exit_code == 0){
        my $cleanup_cmd = "rm -fr $subdir";
        Genome::Sys->shellcmd(cmd => $cleanup_cmd);
        die $self->error_message("Specified git branch ($branch) could not be checked out, repository has been deleted");
    }

    my $latest = $dir . "/" . $source;
    $latest .= "/" . $subsource if $subsource;
    $latest .= "/" . $species if $species;
    $latest .= "/" . "/latest";
    if (-e $latest) {
        unlink $latest;
    }

    symlink "$dir/$subdir", $latest;

    return 1;
}

1;
