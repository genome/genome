package Genome::Db::Command::InstallFromGitRepo;
use strict;
use warnings;
use Genome;

class Genome::Db::Command::InstallFromGitRepo {
    is => 'Command::V2',
    has => [
        branch => {    is => 'Text',
                       doc => 'the branch of database (git repo) to install, corresponding to a branch name in the genome-db-$SOURCE-data repository in genome-vendor on github' },
        source  => {   is => 'Text',
                       is_optional => 1,
                       doc => 'the data source' },
        subsource => { is => 'Text',
                       is_optional => 1,
                       doc => 'the data sub-source (if applicable. e.g. tgi)' },
        species => {   is => 'Text',
                       is_optional => 1,
                       doc => 'the species of the source (if applicable, e.g. human or mouse)' },
        repo    => {   is => 'Text',
                       is_optional => 1,
                       doc => 'override the repository URL' },
    ],
    doc => 'install some versioned data set on the local GMS for which there is a branch in the canonical repository'
};

sub help_detail {
    return <<EOS
Install a file-based database from one of the canonical sets at http://github.com/genome-vendor/genome-db-\$SOURCE-data.git.

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
    my $build = '';

    #If there is 'human-' or 'mouse-' a the beginning of the branch name, do not use this in the directory name (arbitrary convention)
    #Other arbitrary parsing of the final dirname should go here
    my $dirname;
    if(($branch =~ /\w+\-(build\d+)\-(\d+)/) && ($source eq 'dbsnp')){
        $build = $1;
        $dirname = $2;
    }elsif ($branch =~ /human\-(.*)/){
        $dirname = $1;
    }elsif($branch =~ /mouse\-(.*)/){
        $dirname = $1; 
    }else{
        $dirname = $branch;
    }
    #Override 'dbsnp' source name 'dbsnp' to 'genome-db-dbsnp'. We have *different* DBs name 'dbsnp' and 'genome-db-dbsnp'...
    $source = 'genome-db-dbsnp' if ($source eq 'dbsnp');

    $self->status_message("Repo is: " . $repo);
    $self->status_message("Source is: " . $source);
    $self->status_message("Subsource is: " . $subsource) if $self->subsource;
    $self->status_message("Species is: " . $species) if $self->species;
    $self->status_message("Build is: " . $build) if $build;
    $self->status_message("Dirname is: " . $dirname);
    $self->status_message("Branch is: " . $branch);

    my $subdir = $dirname;
    $subdir = $build . "/" . $subdir if $build;
    $subdir = $species . "/" . $subdir if $species;
    $subdir = $subsource . "/" . $subdir if $subsource;
    $subdir = $source . "/" . $subdir;

    my $fulldir = $dir . "/" . $subdir;

    #If the repo/branch directory already exists, do nothing
    if (-e $fulldir && -d $fulldir){
      die $self->error_message("It appears that this branch has already been cloned to $fulldir.  Aborting.  Delete it if you want to clone from scratch.");
    }

    #Attempt to clone the specified git branch into the desired sub-directory of the database directory specified by $ENV{$GENOME_DB}
    my $clone_cmd = "cd $dir; git clone $repo -b $branch $subdir";
    Genome::Sys->shellcmd(cmd => $clone_cmd);

    unless (-e $fulldir && -d $fulldir){
      die $self->error_message("git clone command failed to create the expected directory");
    }

    #Make sure the specified branch exists and was checked out successfully
    my $exit_code = `bash -c 'cd $fulldir; git ls-remote --exit-code . origin/$branch &> /dev/null; echo \$\?'`;
    chomp($exit_code);
    unless ($exit_code == 0){
        my $cleanup_cmd = "rm -fr $fulldir";
        Genome::Sys->shellcmd(cmd => $cleanup_cmd);
        die $self->error_message("Specified git branch ($branch) could not be checked out, repository has been deleted");
    }

    #Create a 'latest' symlink to the db just installed
    my $latest = $dir . "/" . $source;
    $latest .= "/" . $subsource if $subsource;
    $latest .= "/" . $species if $species;
    $latest .= "/" . $build if $build;
    $latest .= "/" . "/latest";
    if (-e $latest) {
        unlink $latest;
    }
    symlink "$dir/$subdir", $latest;

    #If the installed repo contains a 'join_db.pl' there may be large files that were split to go on github that now need to be joined
    my $join_script = $fulldir . "/" . "join_db.pl";
    if (-e $join_script){
      Genome::Sys->shellcmd(cmd => $join_script);
    }

    return 1;
}

1;
