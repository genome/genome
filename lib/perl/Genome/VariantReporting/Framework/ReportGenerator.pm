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
};

sub plan {
    my $self = shift;

    return Genome::VariantReporting::Framework::Plan::MasterPlan->create_from_json($self->plan_json);
}
Memoize::memoize('plan');

sub execute {
    my $self = shift;

    my @entry_processors;
    for my $reporter_plan ($self->plan->reporter_plans) {
        $reporter_plan->object->initialize($self->output_directory);
        push @entry_processors, Genome::VariantReporting::Framework::EntryProcessor->create(
            reporter_plan => $reporter_plan,
            translations => $self->translations,
        );
    }
    $self->debug_message("Reading from: ".$self->input_vcf."\n");
    my $vcf_reader = Genome::File::Vcf::Reader->new($self->input_vcf);
    while (my $entry = $vcf_reader->next) {
        for my $processor (@entry_processors) {
            $processor->process_entry($entry);
        }
    }

    for my $reporter_plan ($self->plan->reporter_plans) {
        $reporter_plan->object->finalize();
    }
    return 1;
}

1;
