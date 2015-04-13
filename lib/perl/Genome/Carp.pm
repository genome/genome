package Genome::Carp;

use strict;
use warnings;

use Carp qw();
use Exporter qw(import);
use Sub::Install qw(install_sub);

my @carp_methods = qw(carpf croakf confessf cluckf longmessf shortmessf);
our @EXPORT_OK = (@carp_methods, 'dief', 'warnf');

for my $method_name (@carp_methods) {
    my $carp_method_name = substr($method_name, 0, -1);
    my $carp_method = Carp->can($carp_method_name);
    unless ($carp_method) {
        die 'method not found: Carp::' . $carp_method_name;
    }

    install_sub({
        code => sub {
            my ($template, @args) = @_;
            local $Carp::CarpLevel = $Carp::CarpLevel + 1;
            return $carp_method->(sprintf($template, @args));
        },
        into => __PACKAGE__,
        as => $method_name,
    });
}

my %method_name = (
    croak => 'dief',
    carp => 'warnf',
);
for my $carp_method_name (keys %method_name) {
    my $method_name = $method_name{$carp_method_name};
    my $carp_method = Carp->can($carp_method_name);
    install_sub({
        code => sub {
            my ($template, @args) = @_;
            return $carp_method->(sprintf($template, @args));
        },
        into => __PACKAGE__,
        as => $method_name,
    });
}

1;
