package Genome::Annotation::Command::AnnotateBuild;

use strict;
use warnings FATAL => 'all';
use Genome;
use Genome::Annotation::Plan;

class Genome::Annotation::Command::AnnotateBuild {
    is => 'Command::V2',
    has_input => [
        build => {
            is => 'Genome::Model::Build',
            doc => "The build who's variants you wish to anotate (Currently ".
                "only supports Somatic Variation and Somatic Validation)",
        },
        variant_type => {
            is => 'Text',
            is_output => 1,
            valid_values => ['snvs', 'indels'],
            doc => "The type of variants used for annotation",
        },
        output_directory => {
            is => 'Path',
            is_output => 1,
            doc => "The location of the reports to be generated",
        },
        plan_file => {
            is => 'Path',
            doc => 'A plan (yaml) file describing the annotation process',
        },
        log_directory => {
            is => 'Path',
            doc => 'The directory where log files will be written.',
        },
    ],
};

sub execute {
    my $self = shift;

    $self->status_message("Constructing plan from file (%s)", $self->plan_file);
    my $plan = Genome::Annotation::Plan->create_from_file($self->plan_file);
    $self->status_message("Validating plan...");
    $plan->validate();
    $self->status_message("Plan is valid.");

    local $ENV{UR_DUMP_DEBUG_MESSAGES} = 1;
    local $ENV{UR_COMMAND_DUMP_DEBUG_MESSAGES} = 1;
    local $ENV{UR_DUMP_STATUS_MESSAGES} = 1;
    local $ENV{UR_COMMAND_DUMP_STATUS_MESSAGES} = 1;
    local $ENV{WF_USE_FLOW} = 1;

    Genome::Annotation::MasterCommand->execute(
        build => $self->build,
        variant_type => $self->variant_type,
        output_directory => $self->output_directory,
        log_directory => $self->log_directory,
        plan => $plan,
    );

    $self->status_message("Writing plan file to output_directory (%s)",
        $self->output_directory);
    $plan->write_to_file(File::Spec->join($self->output_directory, 'plan.yaml'));

    $self->status_message("Annotation complete, reports are located at (%s).",
        $self->output_directory);
}


1;
