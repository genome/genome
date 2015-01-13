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
    return sort_vaf_headers(
        headers      => [$component->create_sample_specific_field_names([vaf_fields()])],
        sample_names => [$component->sample_names],
    );
}

sub per_library_vaf_headers {
    my $component = shift;
    return sort_vaf_headers(
        headers       => [$component->create_library_specific_field_names([vaf_fields()])],
        sample_names  => [$component->sample_names],
        library_names => [$component->library_names],
    );
}

sub sort_vaf_headers {
    my %params = @_;
    return map { $_->{header} } sort {
        $a->{label_priority} <=> $b->{label_priority} ||
        $a->{sample} cmp $b->{sample} ||
        $a->{library} cmp $b->{library} ||
        $a->{vaf_priority} <=> $b->{vaf_priority}
    } annotate_headers(%params);
}

my %VAF_CONVERSIONS = (
    '_ref_count$' => 1,
    '_var_count$' => 2,
    '_vaf$'       => 3,
);

my %LABEL_CONVERSIONS = (
    '^Normal'    => 1,
    '^Discovery' => 2,
    '^Followup'  => 3,
);

sub annotate_headers {
    my %params        = @_;
    my $headers       = $params{headers};
    my $sample_names  = $params{sample_names};
    my $library_names = $params{library_names};

    my @annotated_headers;
    for my $header (@$headers) {
        my %annotated_header = (
            header => $header,
            library => 0,
            label_priority => 0,
            sample => 0,
            vaf_priority => 0,
        );
        while (my ($regex, $sort_priority) = each %VAF_CONVERSIONS) {
            if ($header =~ m/$regex/) {
                $annotated_header{vaf_priority} = $sort_priority;
            };
        }
        for my $library (@$library_names) {
            if ($header =~ m/$library/) {
                $annotated_header{library} = $library;
            }
        }
        for my $sample (@$sample_names) {
            if ($header =~ m/$sample/) {
                $annotated_header{sample} = $sample;
            }
        }
        while (my ($regex, $sort_priority) = each %LABEL_CONVERSIONS) {
            if ($header =~ m/$regex/) {
                $annotated_header{label_priority} = $sort_priority;
            }
        }
        push @annotated_headers, \%annotated_header;
    }
    return @annotated_headers;
}

1;
