package Genome::VariantReporting::Framework::MergeReports;

use strict;
use warnings;
use Genome;

class Genome::VariantReporting::Framework::MergeReports {
    is => 'Genome::Command::DelegatesToResult',
    has_input => [
        report_results => {
            is => 'Genome::VariantReporting::Framework::Component::Report::MergeCompatible',
            doc => 'The reports you wish to merge. They must all be the same type of report (same columns).',
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
            is => 'Genome::VariantReporting::Framework::Component::Report::MergeCompatible',
            doc => 'Use the header from this report_result in the merged report',
            is_optional => 1,
        },
        separator => {
            is => 'Text',
            is_optional => 1,
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
            doc => "An array of strings of the format <Report ID>|<TAG>",
        },
        user => {
            is => 'Genome::Process',
            is_optional => 1,
            id_by => 'process_id',
        },
        process_id => {
            is => 'Text',
        },
    ],
};

sub result_class {
    return "Genome::VariantReporting::Framework::Component::Report::MergedReport";
}

1;
