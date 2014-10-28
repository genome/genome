package Genome::VariantReporting::Suite::BamReadcount::VafInterpreterHelpers;

use strict;
use warnings;
use Genome;

use Exporter 'import';

our @EXPORT_OK = qw(
    basic_available_fields
    basic_field_descriptions
    many_samples_available_fields
    many_samples_field_descriptions
);

sub basic_available_fields {
    return qw(
        vaf
        ref_count
        var_count
        per_library_var_count
        per_library_ref_count
        per_library_vaf
    );
}

sub basic_field_descriptions {
    my $sample_name = shift;
    return (
        vaf => "Variant allele frequency for sample $sample_name",
        ref_count => "Number of reads supporting the reference for sample $sample_name",
        var_count => "Number of reads supporting variant for sample $sample_name",
        per_library_var_count => "Number of reads supporting variant for each library of sample $sample_name",
        per_library_ref_count => "Number of reads supporting the reference for each library of sample $sample_name",
        per_library_vaf => "Variant allele frequency for each library of sample $sample_name",
    );
}

sub many_samples_available_fields {
    my $sample_names = shift;

    my %available_fields = many_samples_field_descriptions($sample_names);
    return keys %available_fields;
}

sub many_samples_field_descriptions {
    my $sample_names = shift;

    my %field_descriptions;
    for my $sample_name (@$sample_names) {
        my %field_descriptions_for_sample = basic_field_descriptions($sample_name);
        for my $field (basic_available_fields()) {
            my $sample_specific_field = Genome::VariantReporting::Framework::Component::WithManySampleNames->create_sample_specific_field_name(
                $field,
                $sample_name
            );
            $field_descriptions{$sample_specific_field} = $field_descriptions_for_sample{$field};
        }
    }
    return %field_descriptions;
}

1;
