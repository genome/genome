package Genome::VariantReporting::Framework::Test::WithTranslationsReport;

use strict;
use warnings FATAL => 'all';
use Genome;
use Data::Dump qw(pp);

class Genome::VariantReporting::Framework::Test::WithTranslationsReport {
    is => 'Genome::VariantReporting::Framework::Component::Report::SingleFile',
};

sub name {
    return '__translated_test__';
}

sub required_interpreters {
    return qw(__translated_test__);
}

sub report {
    my $self = shift;
    my $interpretations = shift;

    $self->_output_fh->print(pp($interpretations) . "\n");
}


1;
