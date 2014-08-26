package Genome::VariantReporting::Dbsnp::Expert;

use strict;
use warnings FATAL => 'all';
use Genome;

class Genome::VariantReporting::Dbsnp::Expert {
    is => 'Genome::VariantReporting::Framework::Component::Expert',
};

sub name {
    'dbsnp';
}

sub run_class {
    my $self = shift;
    my $factory = Genome::VariantReporting::Framework::Factory->create();
    return $factory->get_class('runners', 'joinx');
}


1;
