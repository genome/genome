package Genome::VariantReporting::Framework::Factory;

use strict;
use warnings FATAL => 'all';
use Genome;
use Carp qw(confess);
use Memoize qw();
use File::Spec;
use File::Basename qw(dirname);
use Cwd qw(abs_path);
use Params::Validate qw(validate_pos);

sub search_path {
    return ['Genome::VariantReporting'];
}

use Module::Pluggable
    require => 1,
    search_path => search_path(),
    only => qr(Expert$),
    except => ['Genome::VariantReporting::Component::Expert'],
    sub_name => 'experts';

use Module::Pluggable
    require => 1,
    search_path => search_path(),
    only => qr(Adaptor$),
    except => ['Genome::VariantReporting::Component::Adaptor'],
    sub_name => 'adaptors';

use Module::Pluggable
    require => 1,
    search_path => search_path(),
    only => qr(Run$),
    sub_name => 'runners';


use Module::Pluggable
    require => 1,
    search_path => search_path(),
    only => qr(Filter$),
    except => ['Genome::VariantReporting::Component::Filter'],
    sub_name => 'filters';

use Module::Pluggable
    require => 1,
    search_path => search_path(),
    only => qr(Interpreter$|Filter$),
    except => [
        'Genome::VariantReporting::Component::Filter',
        'Genome::VariantReporting::Component::Interpreter',
    ],
    sub_name => 'interpreters';

use Module::Pluggable
    require => 1,
    search_path => search_path(),
    only => qr(Reporter$),
    sub_name => 'reporters';

class Genome::VariantReporting::Framework::Factory {};

sub names {
    my ($self, $accessor) = @_;
    return keys %{$self->_load($accessor)};
}

sub get_object {
    my ($self, $accessor, $name, $params, $overrides) = validate_pos(@_, 1, 1, 1, 1, 0);

    my $pkg = $self->get_class($accessor, $name);
    return $pkg->create(resolve_params($params, $overrides));
}

sub resolve_params {
    my ($params, $overrides) = @_;
    my %resolved_params = %$params;
    @resolved_params{keys %$overrides} = values %$overrides;
    return %resolved_params;
}

sub get_class {
    my ($self, $accessor, $name) = validate_pos(@_, 1, 1, 1);

    if (exists $self->_load($accessor)->{$name}) {
        my $pkg = $self->_load($accessor)->{$name};
        return $pkg;
    } else {
        confess sprintf("No $accessor with name ($name) available $accessor are:\n    %s\n",
            join("\n    ", $self->names($accessor)));
    }
}

sub _load {
    my ($self, $accessor) = @_;

    my %plugins;
    for my $plugin ($self->$accessor) {
        next unless $plugin->can('name');
        my $name = $plugin->name;
        if (!exists($plugins{$name})) {
            $plugins{$name} = $plugin;
        } else {
            confess sprintf("Attempted to register two $accessor with name (%s):\n    %s",
                $name, join("\n    ", $plugins{$name}, $plugin));
        }
    }
    return \%plugins;
}

1;
