package Genome::VariantReporting::ReportGenerator;

use strict;
use warnings;
use Genome;
use Genome::File::Vcf::Reader;

class Genome::VariantReporting::ReportGenerator {
    is => 'Command::V2',
    has_input => [
        vcf_file => {
            is => 'File',
        },
        plan => {
            is => 'Genome::VariantReporting::Plan',
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

sub execute {
    my $self = shift;

    my @entry_processors;
    for my $reporter_plan ($self->plan->reporter_plans) {
        $reporter_plan->object->initialize($self->output_directory);
        push @entry_processors, Genome::VariantReporting::EntryProcessor->create(
            reporter_plan => $reporter_plan,
            translations => $self->translations,
        );
    }
    $self->debug_message("Reading from: ".$self->vcf_file."\n");
    my $vcf_reader = Genome::File::Vcf::Reader->new($self->vcf_file);
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
