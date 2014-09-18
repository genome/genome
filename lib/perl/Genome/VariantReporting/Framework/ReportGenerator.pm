package Genome::VariantReporting::Framework::ReportGenerator;

use strict;
use warnings;
use Genome;
use Genome::File::Vcf::Reader;
use Memoize;

class Genome::VariantReporting::Framework::ReportGenerator {
    is => 'Command::V2',
    has_input => [
        input_vcf => {
            is => 'File',
        },
        plan_json => {
            is => 'Genome::VariantReporting::Framework::Plan::MasterPlan',
        },
        output_directory => {
            is => 'Path',
            is_output => 1,
        },
        variant_type => {
            is => 'Text',
            valid_values => ['snvs', 'indels'],
        },
        translations => {
            is => 'HASH',
        },
    ],
    has_param => [
        lsf_resource => {
            value => q{-R 'select[mem>16000] rusage[mem=16000]' -M 16000000},
        },
    ],
};

sub plan {
    my $self = shift;

    return Genome::VariantReporting::Framework::Plan::MasterPlan->create_from_json($self->plan_json);
}
Memoize::memoize('plan');

sub execute {
    my $self = shift;

    $self->status_message("Reading from: ".$self->input_vcf."\n");
    my $vcf_reader = Genome::File::Vcf::Reader->new($self->input_vcf);
    while (my $entry = $vcf_reader->next) {
        for my $processor ($self->entry_processors) {
            $processor->process_entry($entry);
        }
    }

    for my $reporter_plan ($self->plan->reporter_plans) {
        $reporter_plan->object->finalize();
    }
    return 1;
}

sub entry_processors {
    my $self = shift;

    my @entry_processors;
    for my $reporter_plan ($self->plan->reporter_plans) {
        my $reporter = $self->object_with_translations($reporter_plan);

        my @filters = map {$self->object_with_translations($_)} $reporter_plan->filter_plans;
        my @interpreters = map {$self->object_with_translations($_)} $reporter_plan->interpreter_plans;

        for my $filter (@filters) {
            $reporter->add_filter_object($filter)
        }
        for my $interpreter (@interpreters) {
            $reporter->add_interpreter_object($interpreter)
        }

        $reporter->initialize($self->output_directory);

        push @entry_processors, Genome::VariantReporting::Framework::EntryProcessor->create(
            reporter => $reporter,
            filters => \@filters,
            interpreters => \@interpreters,
        );
    }
    return @entry_processors;
}
Memoize::memoize('entry_processors');

sub object_with_translations {
    my ($self, $plan) = @_;

    my $object = $plan->object;
    $object->translate_inputs($self->translations);
    return $object;
}


1;
