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
        my $reporter = $reporter_plan->object;
        my @filters = map {$_->object} $reporter_plan->filter_plans;
        my @interpreters = map {$_->object} $reporter_plan->interpreter_plans;

        $reporter->translate_inputs($self->translations);
        for (@filters) {$_->translate_inputs($self->translations)};
        for (@interpreters) {$_->translate_inputs($self->translations)};

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

1;
