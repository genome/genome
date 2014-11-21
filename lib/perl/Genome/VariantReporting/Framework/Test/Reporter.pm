package Genome::VariantReporting::Framework::Test::Reporter;

use strict;
use warnings FATAL => 'all';
use Genome;
use Data::Dump qw(pp);

class Genome::VariantReporting::Framework::Test::Reporter {
    is => 'Genome::VariantReporting::Framework::Component::Reporter::SingleFile',
};

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
