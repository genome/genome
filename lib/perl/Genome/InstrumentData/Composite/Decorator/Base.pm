package Genome::InstrumentData::Composite::Decorator::Base;

use strict;
use warnings;

use Genome;
use Genome::Utility::Text;

class Genome::InstrumentData::Composite::Decorator::Base {
    is_abstract => 1,
    is => 'UR::Singleton',
};

sub decorate {
    my $class = shift;
    my $operation = shift;
    my $params = shift;

    die $class->error_message('Subclass %s must implement decorate()', $class);
}

sub name {
    my $class = shift;

    my @parts = split('::', $class);
    return Genome::Utility::Text::camel_case_to_string($parts[-1], '-');
}

1;
