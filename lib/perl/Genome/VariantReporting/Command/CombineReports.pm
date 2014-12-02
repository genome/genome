package Genome::VariantReporting::Command::CombineReports;

use strict;
use warnings;
use Genome;
use Genome::VariantReporting::Command::CombineReportsResult;

my $REPORT_PKG = $Genome::VariantReporting::Command::CombineReportsResult::REPORT_PKG;

class Genome::VariantReporting::Command::CombineReports {
    is => 'Genome::Command::DelegatesToResult',
    has_input => [
        report_results => {
            is => $REPORT_PKG,
            doc => 'The reports you wish to combine. They must all be the same type of report (same columns).',
            is_many => 1,
        },
        sort_columns => {
            is_optional => 1,
            is => 'Text',
            is_many => 1,
            doc => 'Column names to sort by. If the reports have no headers, provide column numbers (1-based). If not provided, the report will be unsorted.',
        },
        contains_header => {
            is => 'Boolean',
            doc => 'Set to true if the report contains headers'
        },
        use_header_from => {
            is => $REPORT_PKG,
            doc => 'Use the header from this report_result in the combined report',
            is_optional => 1,
        },
        separator => {
            is => 'Text',
            default => "\t",
            doc => 'Field separator for the reports',
        },
        split_indicators => {
            is => 'Text',
            is_optional => 1,
            is_many => 1,
            doc => 'A regular expression that indicates that columns whose headers match should be split up.  These columns contain key:value pairs.  The key will be appended to the column header and the value will be put in the column.  Only valid if the file has a header.',
        },
        entry_sources => {
            is_many => 'Text',
            is => 'Text',
            is_optional => 1,
            doc => "An array of strings of the format <Reporter ID>|<TAG>",
        },
    ],
};

sub result_class {
    return "Genome::VariantReporting::Command::CombineReportsResult";
}

sub input_hash {
    my $self = shift;

    return (
        report_results => [$self->report_results],
        sort_columns => [$self->sort_columns],
        contains_header => $self->contains_header,
        use_header_from => $self->use_header_from,
        separator => $self->separator,
        split_indicators => [$self->split_indicators],
        entry_sources => [$self->entry_sources],
    );
}

1;
