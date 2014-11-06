package Genome::VariantReporting::Framework::Test::WithTranslationsReporter;

use strict;
use warnings FATAL => 'all';
use Genome;
use Data::Dump qw(pp);

class Genome::VariantReporting::Framework::Test::WithTranslationsReporter {
    is => 'Genome::VariantReporting::Framework::Component::Reporter::SingleFile',
};

sub name {
    return '__translated_test__';
}

sub requires_interpreters {
    return qw(__translated_test__);
}

sub report {
    my $self = shift;
    my $interpretations = shift;

    $self->_output_fh->print(pp($interpretations) . "\n");
}


1;
