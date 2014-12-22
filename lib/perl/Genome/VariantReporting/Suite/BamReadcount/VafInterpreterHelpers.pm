package Genome::VariantReporting::Suite::BamReadcount::VafInterpreterHelpers;

use strict;
use warnings;
use Genome;

use Exporter 'import';

our @EXPORT_OK = qw(
    many_samples_field_descriptions
    per_sample_vaf_headers
    per_library_vaf_headers
    many_libraries_field_descriptions
);

sub vaf_fields {
    return qw(
        vaf
        ref_count
        var_count
    );
}

sub per_sample_field_descriptions {
    my $sample_name = shift;
    return (
        vaf => "Variant allele frequency for sample $sample_name",
        ref_count => "Number of reads supporting the reference for sample $sample_name",
        var_count => "Number of reads supporting variant for sample $sample_name",
    );
}

sub per_library_field_descriptions {
    my $library_name = shift;
    return (
        var_count => "Number of reads supporting variant for library $library_name",
        ref_count => "Number of reads supporting the reference for library $library_name",
        vaf => "Variant allele frequency for library $library_name",
    );
}

sub many_libraries_field_descriptions {
    my $component = shift;

    my %field_descriptions;
    for my $library_name ($component->library_names) {
        my %field_descriptions_for_library = per_library_field_descriptions($library_name);
        for my $field (vaf_fields()) {
            my $library_specific_field = $component->create_library_specific_field_name(
                $field,
                $library_name
            );
            $field_descriptions{$library_specific_field} = $field_descriptions_for_library{$field};
        }
    }
    return %field_descriptions;
}

sub many_samples_field_descriptions {
    my $component = shift;

    my %field_descriptions;
    for my $sample_name ($component->sample_names) {
        my %field_descriptions_for_sample = per_sample_field_descriptions($sample_name);
        for my $field (vaf_fields()) {
            my $sample_specific_field = $component->create_sample_specific_field_name(
                $field,
                $sample_name
            );
            $field_descriptions{$sample_specific_field} = $field_descriptions_for_sample{$field};
        }
    }
    return %field_descriptions;
}
sub per_sample_vaf_headers {
    my $component = shift;
    return $component->create_sample_specific_field_names(
        [vaf_fields()],
    );
}

sub per_library_vaf_headers {
    my $component = shift;
    return $component->create_library_specific_field_names(
        [vaf_fields()],
    );
}

1;
