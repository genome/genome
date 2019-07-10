package Genome::Model::Build::Command::AbandonPriorBuildsWithStatus;

use strict;
use warnings;

use Genome;

class Genome::Model::Build::Command::AbandonPriorBuildsWithStatus {
    is => 'Command::V2',
    has => [
        status => {
            is => 'Text',
            doc => 'Status of builds to abandon prior builds',
            valid_values => [ grep $_ ne 'Abandoned', @{ Genome::Model::Build->__meta__->property(property_name => 'status')->valid_values } ],
        },
        models => {
            is => 'Genome::Model',
            is_many => 1,
            shell_args_position => 1,
            doc => 'models whose prior builds to abandon',
            require_user_verify => 1,
        },
    ],
    doc => 'abandon prior builds for a model',
};

sub help_detail {
    return <<EOHELP
Abandon prior builds for a model with a matching status, leaving only the most recent build.

Suppose a model has these builds (using the format from `genome model build status`):

 Model: example-model.prod-cwl   Model ID: abcdef1234567890abcdef1234567890
 BUILD_ID        STATUS
 abcdef1234567890abcdef1234567891        Succeeded
 abcdef1234567890abcdef123456789a        Succeeded
 abcdef1234567890abcdef1234567899        Failed
 abcdef1234567890abcdef1234567894        Failed
 abcdef1234567890abcdef1234567898        Succeeded

Running with the --status=Succeeded option would yield:

 Model: example-model.prod-cwl   Model ID: abcdef1234567890abcdef1234567890
 BUILD_ID        STATUS
 abcdef1234567890abcdef1234567891        Succeeded
 abcdef1234567890abcdef123456789a        Abandoned
 abcdef1234567890abcdef1234567899        Failed
 abcdef1234567890abcdef1234567894        Failed
 abcdef1234567890abcdef1234567898        Abandoned

Running instead with the --status=Failed option would yield:

 Model: example-model.prod-cwl   Model ID: abcdef1234567890abcdef1234567890
 BUILD_ID        STATUS
 abcdef1234567890abcdef1234567891        Succeeded
 abcdef1234567890abcdef1234567892        Succeeded
 abcdef1234567890abcdef1234567893        Failed
 abcdef1234567890abcdef1234567894        Abandoned
 abcdef1234567890abcdef1234567895        Succeeded

Running twice, once with each of the above, will leave only one build with each status:

 Model: example-model.prod-cwl   Model ID: abcdef1234567890abcdef1234567890
 BUILD_ID        STATUS
 abcdef1234567890abcdef1234567891        Succeeded
 abcdef1234567890abcdef123456789a        Abandoned
 abcdef1234567890abcdef1234567899        Failed
 abcdef1234567890abcdef1234567894        Abandoned
 abcdef1234567890abcdef1234567898        Abandoned

EOHELP
}

sub execute {
    my $self = shift;

    my @builds_to_abandon;

    my $status = $self->status;

    for my $m ($self->models) {
        my @builds = $m->builds;

        my $found = 0;
        for my $build (@builds) {
            next unless $build->status eq $status;

            if ($found) {
                push @builds_to_abandon, $build;
            } else {
                $found = 1;
            }
        }
    }

    unless (@builds_to_abandon) {
        $self->status_message("Found no matching builds to abandon.");
        return 1;
    }

    my $abandon_cmd = Genome::Model::Build::Command::Abandon->create(
        builds => \@builds_to_abandon,
    );

    $abandon_cmd->builds([ $abandon_cmd->_limit_results_for_builds(@builds_to_abandon) ]);

    return $abandon_cmd->execute;
}

1;
