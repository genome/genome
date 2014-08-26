package Genome::VariantReporting::Joinx::Expert;

use strict;
use warnings FATAL => 'all';
use Genome;

class Genome::VariantReporting::Joinx::Expert {
    is => 'Genome::VariantReporting::Framework::Component::Expert',
    is_abstract => 1,
};

sub run_class {
    my $self = shift;
    my $factory = Genome::VariantReporting::Framework::Factory->create();
    return $factory->get_class('runners', 'joinx');
}


1;
