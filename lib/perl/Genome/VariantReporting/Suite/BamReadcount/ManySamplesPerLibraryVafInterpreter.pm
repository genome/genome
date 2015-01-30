package Genome::VariantReporting::Suite::BamReadcount::ManySamplesPerLibraryVafInterpreter;

use strict;
use warnings;
use Genome;
use Genome::VariantReporting::Suite::BamReadcount::VafInterpreter;
use Genome::VariantReporting::Suite::BamReadcount::VafInterpreterHelpers qw(
    many_libraries_field_descriptions
);

class Genome::VariantReporting::Suite::BamReadcount::ManySamplesPerLibraryVafInterpreter {
    is => [
        'Genome::VariantReporting::Framework::Component::Interpreter',
        'Genome::VariantReporting::Framework::Component::WithManySampleNames',
        'Genome::VariantReporting::Framework::Component::WithManyLibraryNames',
    ],
    has => [],
    doc => 'Calculate the variant allele frequency, number of reads supporting the reference, and number of reads supporting variant for the libraries of multiple samples',
};

sub name {
    return 'many-samples-per-library-vaf';
}

sub requires_annotations {
    return ('bam-readcount');
}

sub field_descriptions {
    my $self = shift;
    return many_libraries_field_descriptions($self);
}

sub _interpret_entry {
    my $self = shift;
    my $entry = shift;
    my $passed_alt_alleles = shift;

    my %return_values;
    for my $sample_name ($self->sample_names) {
        my $interpreter = Genome::VariantReporting::Suite::BamReadcount::PerLibraryVafInterpreter->create(
            sample_name => $sample_name,
            sample_name_label => $self->sample_name_labels->{$sample_name},
            library_names => [$self->library_names],
            library_name_labels => $self->library_name_labels,
        );
        my %result = $interpreter->interpret_entry($entry, $passed_alt_alleles);
        for my $alt_allele (@$passed_alt_alleles) {
            for my $field_name (keys %{$result{$alt_allele}}) {
                $return_values{$alt_allele}->{$field_name} = $result{$alt_allele}->{$field_name};
            }
        }
    }
    return %return_values;
}

1;
