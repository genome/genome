package Genome::Model::TestPipeline;

use strict;
use warnings;
use Genome;

class Genome::Model::TestPipeline {
    is => 'Genome::ModelDeprecated',
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
        # NOTE: these are made up parameters just as examples
        # A processing profile shouldn't really have params that specify shell commands :)
        some_command_name => {
            doc => 'the name of a single command to run',
            valid_values => ['ls','cat','wc'],
            default_value => ['ls'],
        },
        some_args => {
            is_optional => 1,
            doc => 'the arguments to use',
        }
    ],
    doc => "an example processing profile which runs one shell command and catches its output"
};

sub _execute_build {
    my ($self, $build) = @_;
    warn "executing build logic for " . $self->__display_name__ . ':' .  $build->__display_name__ . "\n";

    # combine params with build inputs and produce output in the build's data directory

    my $cmd = $self->some_command_name;
    my $args = $self->some_args;

    my @inputs = $build->inputs();

    my $dir = $build->data_directory;

    my $exit_code = system "$cmd $args @inputs >$dir/output 2>$dir/errors";
    $exit_code = $exit_code >> 8;
    if ($exit_code != 0) {
        $self->error_message("Failed to run $cmd with args $args!  Exit code: $exit_code.");
        return;
    }

    return 1;
}

sub _validate_build {
    my $self = shift;
    my $dir = $self->data_directory;

    my @errors;
    unless (-e "$dir/output") {
        my $e = $self->error_message("No output file $dir/output found!");
        push @errors, $e;
    }
    unless (-e "$dir/errors") {
        my $e = $self->error_message("No output file $dir/errors found!");
        push @errors, $e;
    }

    if (@errors) {
        return;
    }
    else {
        return 1;
    }
}

1;

