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
        config => { is => 'Text', is_optional => 1 },
    ],
};

sub extend_class_with_config {
    my ($class, $desc) = @_;

    while (my ($prop_name, $prop_desc) = each(%{ $desc->{has} })) {
        next unless (exists $prop_desc->{config} && $prop_desc->{config});

        my $spec = Genome::Config::spec($prop_desc->{config});

        my $key = $spec->key;
        my $default_method_name = '__default_' . $prop_name . '__';
        my $default = sub { Genome::Config::get($key) };
        install_sub({
            code => $default,
            into => $prop_desc->{class_name},
            as   => $default_method_name,
        });

        set_prop_desc_with_warn($prop_desc, 'data_type', $spec->type);
        set_prop_desc_with_warn($prop_desc, 'calculated_default', $default_method_name);
    }

    return $desc;
}

sub set_prop_desc_with_warn {
    my ($prop_desc, $key, $value) = @_;
    if ($prop_desc->{$key}) {
        warnf '%s has %s but overriding from config', $prop_desc->{class_name}, $key;
    }
    $prop_desc->{$key} = $value;
}

1;
