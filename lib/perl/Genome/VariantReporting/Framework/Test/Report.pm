package Genome::VariantReporting::Framework::Test::Report;

use strict;
use warnings FATAL => 'all';
use Genome;
use Data::Dump qw(pp);

class Genome::VariantReporting::Framework::Test::Report {
    is => 'Genome::VariantReporting::Framework::Component::Report::SingleFile',
    has_optional_transient => [
        _report_path => {
            is_structural => 1,
            is => 'Text',
        }
    ],
};

sub report_path {
    my $self = shift;
    if (defined($self->_report_path)) {
        return $self->_report_path;
    } else {
        return $self->SUPER::report_path(@_);
    }
}

sub name {
    return '__test__';
}

sub required_interpreters {
    return qw(__test__);
}

sub report {
    my $self = shift;
    my $interpretations = shift;

    $self->_output_fh->print(pp($interpretations) . "\n");
}


1;
