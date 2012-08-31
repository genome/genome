package Genome::Sys::Base;
use strict;
use warnings;
use Genome;

class Genome::Sys::Base {
    doc => 'a base class for classes which want Genome::Sys functionality as methods',
};

sub CAN {
    my ($class_name, $method, $self, @params) = @_;
    my $m = Genome::Sys->can($method);
    return unless $m;
    my $full_name = join( '::', $class_name, $method); 
    my $code =  sub { 
                    shift; 
                    return Genome::Sys->$method(@_) 
                };
    Sub::Install::reinstall_sub({
        into => $class_name,
        as   => $method,
        code => Sub::Name::subname($full_name => $code) 
    });

    return $code;
}

# this must be done at run time so the class is defined first
eval "use Class::AutoloadCAN";

1;

