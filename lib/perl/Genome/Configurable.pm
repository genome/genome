package Genome::Configurable;

use strict;
use warnings;

use Genome::Carp qw(warnf);
use Genome::Config qw();
use Sub::Install qw(install_sub);
use UR qw();

class Genome::Configurable {
    is_abstract => 1,
    subclass_description_preprocessor => __PACKAGE__ . '::extend_class_with_config',
    attributes_have => [
        config => { is => 'Text', is_optional => 1, is_many => 1 },
    ],
};

sub extend_class_with_config {
    my ($class, $desc) = @_;

    while (my ($prop_name, $prop_desc) = each(%{ $desc->{has} })) {
        next unless (exists $prop_desc->{config} && $prop_desc->{config});

        my $key = $prop_desc->{config};
        my $default;
        if (ref($key) && ref($key) eq 'ARRAY') {
            my @keys = @{$key};
            my @specs = map { Genome::Config::spec($_) } @keys;
            $default = sub { Genome::Config::get_first(@keys) };
        }
        else {
            my $spec = Genome::Config::spec($key);
            $default = sub { Genome::Config::get($key) };
        }

        my $default_method_name = '__default_' . $prop_name . '__';
        install_sub({
            code => $default,
            into => $prop_desc->{class_name},
            as   => $default_method_name,
        });

        set_prop_desc_with_warn($prop_desc, 'calculated_default', $default_method_name);
    }

    return $desc;
}

sub set_prop_desc_with_warn {
    my ($prop_desc, $key, $value) = @_;
    if (exists $prop_desc->{$key}) {
        warnf '%s has %s but overriding from config', $prop_desc->{class_name}, $key;
    }
    $prop_desc->{$key} = $value;
}

1;

=pod

=head1 NAME

Genome::Configurable

=head1 DESCRIPTION

Genome::Configurable is a role that injects the C<calculated_default> attribute
into UR class properties based on the C<config> attribute such that the default
value is determined via L<Genome::Config> when the object is created.

Genome::Configurable is primarily intended for commands to incorporate command
line arguments as a new primary source above the L<Genome::Config> API.

=cut
