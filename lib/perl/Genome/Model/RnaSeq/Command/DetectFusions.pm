package Genome::Model::RnaSeq::Command::DetectFusions;

use strict;
use warnings;

use Genome;

class Genome::Model::RnaSeq::Command::DetectFusions {
    is => "Genome::Command::Base",
    has_input => [
        detector_name => {
            is => 'Text',
            doc => 'the name of the fusion detector to run',
        },
        detector_version => {
            is => 'Text',
            doc => 'the version of the fusion detector to run',
        },
        detector_params => {
            is => 'Text',
            doc => 'parameters for the chosen fusion detector',
        },
        build_id => {
            is => 'Text',
        },
        build => {
            is => "Genome::Model::Build::RnaSeq",
            id_by => 'build_id',
        },
    ],
    doc => 'run a specified fusion detector',
};

# Keys: detector_name
my %COMMANDS = (
    'chimerascan' => 'Genome::Model::RnaSeq::Command::DetectFusions::Chimerascan::FixedReadLength',
    'chimerascan-vrl' => 'Genome::Model::RnaSeq::Command::DetectFusions::Chimerascan::VariableReadLength',
);

sub execute {
    my $self = shift;

    my $command_class = $COMMANDS{$self->detector_name};
    my $cmd = $command_class->create(
        detector_version => $self->detector_version,
        detector_params => $self->detector_params,
        build => $self->build,
    );

    my $rv = $cmd->execute();
    if ($cmd->result) {
        $self->_link_build_to_result($cmd->result);
    }

    return $rv;
}

sub _link_build_to_result {
    my $self = shift;
    my $result = shift;

    Genome::Sys->create_symlink($result->output_dir, $self->build->data_directory . "/fusions");
    my $link = $result->add_user(user => $self->build, label => 'uses');
    if ($link) {
        $self->status_message("Linked result " . $result->id . " to the build");
    }
    else {
        $self->error_message(
            "Failed to link the build to the result "
            . $result->__display_name__ . "!"
        );
        die $self->error_message;
    }
    return 1;
}

1;
