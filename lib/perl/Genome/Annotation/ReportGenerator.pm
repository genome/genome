package Genome::Annotation::ReportGenerator;

use strict;
use warnings;
use Genome;
use Genome::File::Vcf::Reader;

class Genome::Annotation::ReportGenerator {
    is => 'Command::V2',
    has_input => [
        vcf_file => {
            is => 'File',
        },
        plan => {
            is => 'Genome::Annotation::Plan',
        },
        output_directory => {
            is => 'Path',
            is_output => 1,
        },
        variant_type => {
            is => 'Text',
            valid_values => ['snv', 'indel'],
        },
    ],
};

sub execute {
    my $self = shift;
    for my $reporter_plan ($self->plan->reporter_plans) {
        $reporter_plan->object->initialize($self->output_directory);
    }

    my $vcf_reader = Genome::File::Vcf::Reader->new($self->vcf_file);
    while (my $entry = $vcf_reader->next) {
        for my $reporter_plan ($self->plan->reporter_plans) {
            $reporter_plan->process_entry($entry);
        }
    }
    for my $reporter_plan ($self->plan->reporter_plans) {
        $reporter_plan->object->finalize();
    }
    return 1;
}

1;
