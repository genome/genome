use strict;
use warnings;
use Genome;

class Genome::Model::ImportedAssembly {
    is => 'Genome::ModelDeprecated',
    doc => 'imported assembly',
    has => [
        server_dispatch => {
            is_constant => 1,
            is_class_wide => 1,
            value => 'inline',
            doc => 'lsf queue to submit the launcher or \'inline\''
        },
        job_dispatch => {
            is_constant => 1,
            is_class_wide => 1,
            value => 'inline',
            doc => 'lsf queue to submit jobs or \'inline\' to run them in the launcher'
        }
    ],
    has_param => [
        assembler_name => {
            doc => 'Assembler used to generate the assembly',
        },
        assembler_version => {
            is_optional => 1,
            doc => 'Version of the assembler used to generate the assembly',
        },
        assembly_description => {
            is_optional => 1,
            doc => 'Special description of the assembly',
        },
    ],
};

sub _execute_build {
    my ($self, $build) = @_;

    unless (-d $build->data_directory) {
    $self->error_message("Failed to find build directory: ".$build->data_directory);
    return;
    }
    else {
    $self->status_message("Created build directory: ".$build->data_directory);
    }

    return 1;
}

1;

