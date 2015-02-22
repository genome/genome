package Genome::VariantReporting::Command::CreateReport;

use strict;
use warnings FATAL => 'all';
use Genome::VariantReporting::Framework::Dag qw(generate_dag);
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
        plan_file => {
            is => 'Path',
            doc => 'A plan (yaml) file describing the report generation workflow',
        },
        translations_file => {
            is => 'Path',
            doc => 'A yaml file containing key-value pairs where the key is a value from the plan file that needs to be translated at runtime',
        },
    ],
    has_transient_optional => [
        plan => {
            is => 'Genome::VariantReporting::Framework::Plan::MasterPlan',
        },
        provider => {
            is => 'Genome::VariantReporting::Framework::Component::RuntimeTranslations',
        },
        dag => {
            is => 'Genome::WorkflowBuilder::DAG',
        }
    ],
};

sub process_class {
    return "Genome::VariantReporting::Process::CreateReport";
}

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

    my $p = $self->process_class->create(
        input_vcf => $self->input_vcf,
        variant_type => $self->variant_type,
        report_names => [map {$_->name} $self->plan->report_plans],
    );
    $p->save_plan_file($self->plan_file);
    $p->save_translations_file($self->translations_file);

    $self->status_message("Constructing workflow from plan file (%s)", $self->plan_file);
    $p->run(workflow_xml => $self->dag->get_xml,
        workflow_inputs => $self->workflow_inputs($p->id),
    );

    return $p;
}

sub workflow_inputs {
    my $self = shift;
    my $process_id = shift;
    return {
        process_id => $process_id,
        %{$self->dag->constant_values},
    };
}

sub plan {
    my $self = shift;

    unless (defined($self->__plan)) {
        $self->debug_message("Constructing plan from file (%s)", $self->plan_file);
        my $plan = Genome::VariantReporting::Framework::Plan::MasterPlan->create_from_file($self->plan_file);
        $self->debug_message("Validating plan...");
        $plan->validate();
        $self->debug_message("Plan is valid.");

        $self->debug_message("Checking for compatibility between translations and plan...");
        $plan->validate_translation_provider($self->provider);
        $self->debug_message("Translations file is compatible with plan.");

        $self->debug_message("Translating plan");
        $plan->translate($self->provider->translations);

        $self->__plan($plan);
    }
    return $self->__plan;
}

sub provider {
    my $self = shift;

    unless (defined($self->__provider)) {
        $self->debug_message("Constructing translation-provider from file (%s)", $self->translations_file);
        my $provider = Genome::VariantReporting::Framework::Component::RuntimeTranslations->create_from_file($self->translations_file);

        $self->__provider($provider);
    }
    return $self->__provider;
}

sub dag {
    my $self = shift;

    unless (defined($self->__dag)) {
        my $dag = generate_dag($self->plan, $self->variant_type);
        $dag->declare_constant(
            input_vcf => $self->input_vcf,
            variant_type => $self->variant_type,
            plan_json => $self->plan->as_json,
        );
        $self->__dag($dag);
    }
    return $self->__dag;
}

1;
