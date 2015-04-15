package Genome::InstrumentData::Composite::Decorator;

use strict;
use warnings;

use Genome;
use Genome::Utility::Text;

class Genome::InstrumentData::Composite::Decorator {
    is => 'UR::Singleton',
};

sub decorate {
    my $class = shift;
    my $operation = shift;
    my $decoration = shift;

    my $decorator = $class->_resolve_decorator($decoration);
    return $decorator->decorate($operation, $decoration->{params});
}

sub _resolve_decorator {
    my $class = shift;
    my $decoration = shift;

    my $name = $decoration->{name};
    my $decorator_class = join('::', __PACKAGE__, Genome::Utility::Text::string_to_camel_case($name, '-'));

    my $type = UR::Object::Type->get(class_name => $decorator_class);
    unless($type) {
        die $class->error_message('No decorator found for %s.', $name);
    }

    return $decorator_class;
}

1;
