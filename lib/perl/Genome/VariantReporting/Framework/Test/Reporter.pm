package Genome::VariantReporting::Framework::Test::Reporter;

use strict;
use warnings FATAL => 'all';
use Genome;
use Data::Dump qw(pp);

class Genome::VariantReporting::Framework::Test::Reporter {
    is => 'Genome::VariantReporting::Framework::Component::Reporter::SingleFile',
    has_optional_transient => [
        _report_path => {
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
