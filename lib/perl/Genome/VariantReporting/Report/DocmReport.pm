package Genome::VariantReporting::Report::DocmReport;

use strict;
use warnings;
use Genome;
use Genome::VariantReporting::Suite::BamReadcount::VafInterpreterHelpers qw(
    per_sample_vaf_headers
    per_library_vaf_headers
    );

class Genome::VariantReporting::Report::DocmReport {
    is => [ 'Genome::VariantReporting::Report::WithHeader',
        'Genome::VariantReporting::Framework::Component::WithManyLibraryNames',
        'Genome::VariantReporting::Framework::Component::WithManySampleNames'],
    has_input => [
        sample_names => {
            is => 'Text',
            is_many => 1,
            is_translated => 1,
            doc => 'List of sample names to be used in the report',
        },
        library_names => {
            is => 'Text',
            is_many => 1,
            is_translated => 1,
            doc => 'List of library names to be used in the report',
        },
    ],
    has_transient_optional_translated => [
        sample_name_labels => {
            is => 'HASH',
            default => {},
            doc => 'Hash of sample_name to label',
        },
        library_name_labels => {
            is => 'HASH',
            default => {},
            doc => 'Hash of library_name to label',
        },
    ],
    doc => 'Output readcount information from bam readcount',
};

sub name {
    return 'docm';
}

sub required_interpreters {
    return qw(position vaf per-library-vaf);
}

sub headers {
    my $self = shift;
    my @headers = qw/
        chromosome_name
        start
        stop
        reference
        variant
    /;
    push @headers, per_sample_vaf_headers($self);
    push @headers, per_library_vaf_headers($self);

    return @headers;
}

1;
