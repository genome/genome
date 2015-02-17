package Genome::VariantReporting::Framework::MergeReports;

use strict;
use warnings;
use Genome;

class Genome::VariantReporting::Framework::MergeReports {
    is => 'Genome::Command::DelegatesToResult',
    has_input => [
        base_report => {
            is => 'Genome::VariantReporting::Framework::Component::Report::MergeCompatible',
            doc => 'The main report to be merged.'
        },
        supplemental_report => {
            is => 'Genome::VariantReporting::Framework::Component::Report::MergeCompatible',
            doc => 'The report that is to be added to the main report.'
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
        separator => {
            is => 'Text',
            is_optional => 1,
            default => "\t",
            doc => 'Field separator for the reports',
        },
        entry_sources => {
            is_many => 'Text',
            is => 'Text',
            is_optional => 1,
            doc => "An array of strings of the format <Report ID>|<TAG>",
        },
        process_id => {
            is => 'Text',
        },
    ],
    has_transient_optional => [
        requestor => {
            is => 'Genome::Process',
            id_by => 'process_id',
        },
        user => {
            is => 'Genome::Process',
            id_by => 'process_id',
        },
    ],
};

sub result_class {
    return "Genome::VariantReporting::Framework::Component::Report::MergedReport";
}

1;
