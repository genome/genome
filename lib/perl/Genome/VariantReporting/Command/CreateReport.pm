package Genome::VariantReporting::Command::CreateReport;

use strict;
use warnings FATAL => 'all';
use Genome::VariantReporting::Framework::Dag qw(generate_dag);
use Memoize qw();
use Genome;

class Genome::VariantReporting::Command::CreateReport {
    is => 'Command::V2',
    has_input => [
        input_vcf => {
            is => 'Path',
            doc => 'The source of variants to create a report from',
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
            doc => 'A plan (yaml) file describing the report generation workflow',
        },
        translations_file => {
            is => 'Path',
            doc => 'A yaml file containing key-value pairs where the key is a value from the plan file that needs to be translated at runtime',
        },
        log_directory => {
            is => 'Path',
            doc => 'The directory where log files will be written.',
        },
    ],
};

sub __errors__ {
    my $self = shift;
    my @errors = $self->SUPER::__errors__(@_);

    unless (-s $self->input_vcf) {
        push @errors, UR::Object::Tag->create(
            type => 'invalid',
            properties => [qw/ input_vcf /],
            desc => sprintf('Input vcf specified (%s) does not exist!', $self->input_vcf),
        );
    }

    return @errors;
}

sub execute {
    my $self = shift;

    $self->status_message("Executing workflow.");
    $self->dag->execute(
        $self->params_for_execute
    );

    $self->status_message("Writing plan file and provider file to output_directory (%s)",
        $self->output_directory);
    $self->plan->write_to_file(File::Spec->join($self->output_directory, 'plan.yaml'));
    $self->provider->write_to_file(File::Spec->join($self->output_directory, 'resources.yaml'));

    $self->status_message("Report Generation complete, reports are located at (%s).",
        $self->output_directory);

    return 1;
}

sub params_for_execute {
    my $self = shift;
    return (
        input_vcf => $self->input_vcf,
        variant_type => $self->variant_type,
        output_directory => $self->output_directory,
        plan_json => $self->plan->as_json,
        provider_json => $self->provider->as_json,
        translations => $self->provider->translations,
    );
}

sub plan {
    my $self = shift;

    $self->status_message("Constructing plan from file (%s)", $self->plan_file);
    my $plan = Genome::VariantReporting::Framework::Plan::MasterPlan->create_from_file($self->plan_file);
    $self->status_message("Validating plan...");
    $plan->validate();
    $self->status_message("Plan is valid.");
    return $plan;
}
Memoize::memoize('plan');

sub provider {
    my $self = shift;

    $self->status_message("Constructing translation-provider from file (%s)", $self->translations_file);
    my $provider = Genome::VariantReporting::Framework::Component::RuntimeTranslations->create_from_file($self->translations_file);

    $self->status_message("Checking for compatibility between translations and plan...");
    $self->plan->validate_translation_provider($provider);
    $self->status_message("Translations file is compatible with plan.");
    return $provider;
}
Memoize::memoize('provider');

sub dag {
    my $self = shift;

    $self->status_message("Constructing workflow from plan.");
    my $dag = generate_dag($self->plan, $self->variant_type);

    $self->status_message("Setting log-directory to (%s)", $self->log_directory);
    Genome::Sys->create_directory($self->log_directory);
    $dag->log_dir($self->log_directory);

    return $dag;
}
Memoize::memoize('dag');

1;
