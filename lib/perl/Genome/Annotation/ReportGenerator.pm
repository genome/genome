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
            my $filter_results = initialize_filters($entry);
            for my $filter_plan ($reporter_plan->filter_plans) {
                combine($filter_results, {$filter_plan->object->process_entry($entry)});
            }
            my @passed_alleles = grep {$filter_results->{$_} == 1} keys %$filter_results;
            my %interpretations;
            for my $interpreter_plan ($reporter_plan->interpreter_plans) {
                $interpretations{$interpreter_plan->object->name} = {$interpreter_plan->object->process_entry($entry, \@passed_alleles)};
            }
            $reporter_plan->object->report(\%interpretations);
        }
    }
    for my $reporter_plan ($self->plan->reporter_plans) {
        $reporter_plan->object->finalize();
    }
    return 1;
}

sub initialize_filters {
    my $entry = shift;
    my %filter_values;
    for my $allele (@{$entry->{alternate_alleles}}) {
        $filter_values{$allele} = 1;
    }
    return \%filter_values;
}

sub combine {
    my $accumulator = shift;
    my $new_result = shift;
    for my $allele (keys %$accumulator) {
        $accumulator->{$allele} = $accumulator->{$allele} & $new_result->{$allele};
    }
    return $accumulator;
}

1;

