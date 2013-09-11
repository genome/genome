package Genome::Model::Build::Command::Validate;

use strict;
use warnings;
use Genome;

class Genome::Model::Build::Command::Validate {
    is => 'Genome::Command::Base',
    doc => "Validate a build.",
    has_input => [
        builds => {
            is_output => 1,
            shell_args_position => 1,
            is_many => 1,
            is => 'Genome::Model::Build',
            doc => 'the build(s) to validate',
        },
        builds_can => {
            is => 'String',
            is_many => 1,
            is_optional => 1,
            doc => 'If set, make sure the builds all "can" do these things ($build->can("method"))' ,
        },
    ],
};

sub execute {
    my $self = shift;

    my $failed = 0;
    if ($self->builds_can) {
        $failed += $self->_check_builds_can;
    }

    if ($failed) {
        die $self->error_message("One or more builds failed validation");
    } else {
        $self->status_message("All builds succeeded validation");
    }

    return 1;
}

sub _check_builds_can {
    my $self = shift;

    my $failed = 0;
    for my $build ($self->builds) {
        for my $method ($self->builds_can) {
            unless ($build->can($method)) {
                $self->error_message($build->id . " is unable to $method");
                $failed = 1;
            }
        }
    }

    return $failed;
}

1;

