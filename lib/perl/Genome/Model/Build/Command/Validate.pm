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
        status => {
            is => 'String',
            is_many => 1,
            is_optional => 1,
            valid_values => [qw(Succeeded Failed Scheduled Abandoned Running Unstartable)],
            doc => 'If set, builds must have one of the supplied statuses',
        },
        reference_sequence => {
            is => 'Genome::Model::Build::ReferenceSequence',
            is_many => 1,
            is_optional => 1,
            doc => 'If set, builds must have one of these reference sequence builds',
        },
    ],
};

sub execute {
    my $self = shift;

    my $failed = 0;
    if ($self->builds_can) {
        $failed += $self->_check_builds_can;
    }
    if ($self->status) {
        $failed += $self->_check_builds_status;
    }
    if ($self->reference_sequence) {
        $failed += $self->_check_builds_reference_sequence;
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

sub _check_builds_status {
    my $self = shift;

    my $failed = 0;
    my @statuses = $self->status;
    for my $build ($self->builds) {
        unless (grep {$_ eq $build->status} @statuses) {
            $self->error_message(sprintf("Status %s not valid for build: %s (must be one of %s)",
                $build->status, $build->id, join(', ', @statuses)));
            $failed = 1;
        }
    }

    return $failed;
}

sub _check_builds_reference_sequence {
    my $self = shift;

    my $failed = 0;
    my @reference_sequences = $self->reference_sequence;
    for my $build ($self->builds) {
        my $rfb;
        if ($build->can('reference_sequence_build')) {
            $rfb = $build->reference_sequence_build;
        } else {
            $self->error_message(sprintf("Build %s has no property 'reference_sequence_build'", $build->id));
            $failed = 1;
            next;
        }

        unless (grep {$_->id eq $rfb->id} @reference_sequences) {
            $self->error_message(sprintf("Reference sequence: %s not valid for build: %s (must be one of %s)",
                $rfb->name, $build->id, join(', ', map {$_->__display_name__} @reference_sequences)));
            $failed = 1;
        }
    }

    return $failed;
}

1;

