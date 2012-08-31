
package Genome::Model::Tools::Galaxy::Setup;

use strict;
use warnings;
use Genome;
use File::Copy;

class Genome::Model::Tools::Galaxy::Setup {
    is  => 'Command::V2',
    has => [
        path => {
            is  => 'Text',
            is_optional => 1,
            shell_args_position => 1,
            doc => 'Galaxy setup path.  Defaults to "galaxy" in your home directory.',
        }
    ],
    doc => 'setup the galaxy software on your system',
};

sub sub_command_sort_position { 1 }

sub execute {
    my $self = shift;

    my $path = $self->path;
    if (!defined($path)) {
        $path = $ENV{HOME} . "/galaxy/";
    }
    my $command = "hg clone https://bitbucket.org/galaxy/galaxy-dist $path";
    $self->status_message("Cloning galaxy from remote repository. This is a 200MB download and may take several minutes");
    system($command);
    unless ($? == 0) {
        $self->warning_message("Encountered non zero exit. Error encountered in cloning Galaxy");
        die();
    }

    $self->status_message("Galaxy has been copied to $path. Installing Genome commands.");
    copy($path . "/tool_conf.xml.sample", $path . "/tool_conf.xml");

    $self->status_message("Updating to the latest revision...");
    my $result = Genome::Model::Tools::Galaxy::Update->execute(
        path => $path,
        pull => 0
    );
    unless ($result) {
        $self->error_message("error updating Galaxy!: " . Genome::Model::Tools::Galaxy::Update->error_message());
    }

    $self->status_message("Start galaxy by running: $path/run.sh...");
    return 1;
}

