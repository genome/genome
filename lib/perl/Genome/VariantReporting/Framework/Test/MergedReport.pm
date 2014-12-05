package Genome::VariantReporting::Framework::Test::MergedReport;

use strict;
use warnings FATAL => 'all';
use Genome;

class Genome::VariantReporting::Framework::Test::MergedReport {
    is => 'Genome::VariantReporting::Framework::Component::Report::MergedReport',
    has_optional_transient => [
        _report_path => {
            is => 'Text',
        }
    ],
};

sub report_path {
    my $self = shift;
    return $self->_report_path;
}


1;
