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
    except => ['Genome::VariantReporting::Framework::Component::Expert',
               'Genome::VariantReporting::Suite::Joinx::Expert'],
    sub_name => 'experts';

use Module::Pluggable
    require => 1,
    search_path => search_path(),
    only => qr(Adaptor$),
    except => ['Genome::VariantReporting::Framework::Component::Adaptor',
               'Genome::VariantReporting::Suite::Joinx::Adaptor'],
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
    except => ['Genome::VariantReporting::Framework::Component::Filter'],
    sub_name => 'filters';

use Module::Pluggable
    require => 1,
    search_path => search_path(),
    only => qr(Interpreter$|Filter$),
    except => [
        'Genome::VariantReporting::Framework::Component::Filter',
        'Genome::VariantReporting::Framework::Component::Interpreter',
    ],
    sub_name => 'interpreters';

use Module::Pluggable
    require => 1,
    search_path => search_path(),
    only => qr(Report$),
    except => [
        'Genome::VariantReporting::Framework::Component::Report',
        'Genome::VariantReporting::Framework::Component::Report::SingleFile',
    ],
    sub_name => 'reports';

class Genome::VariantReporting::Framework::Factory {};

sub names {
    my ($self, $accessor) = @_;
    return keys %{$self->_load($accessor)};
}

sub get_object {
    my ($self, $accessor, $name, $params) = validate_pos(@_, 1, 1, 1, 1);
    my $pkg = $self->get_class($accessor, $name);
    return $pkg->create(%$params);
}

{
    package Genome::VariantReporting::Framework::Factory::Dummy;
    require Carp;

    sub AUTOLOAD {
        my $self = shift;

        my $target_sub = our $AUTOLOAD;
        $target_sub =~ s/.+:://;
        if(exists( $self->{$target_sub} )) {
            return $self->{$target_sub};
        } else {
            my $pkg = $self->{class};
            my $sub = $pkg->can($target_sub);
            Carp::croak "Subroutine $target_sub not found on package $pkg" unless $sub;

            $sub->($self, @_);
        }
    }
}

sub get_dummy_object {
    my ($self, $accessor, $name) = validate_pos(@_, 1, 1, 1);
    my $pkg = $self->get_class($accessor, $name);
    my $params = $pkg->_get_dummy_params;

    $params->{class} = $pkg;
    bless $params, 'Genome::VariantReporting::Framework::Factory::Dummy';

    return $params;
}

sub get_class {
    my ($self, $accessor, $name) = validate_pos(@_, 1, 1, 1);

    my $type_name = _truncate_name($name);

    my $pkg = $self->_load($accessor)->{$type_name};
    if (defined($pkg)) {
        return $pkg;
    } else {
        confess sprintf("No $accessor with name ($type_name) available $accessor are:\n    %s\n",
            join("\n    ", $self->names($accessor)));
    }
}
Memoize::memoize('get_class', LIST_CACHE => 'MERGE');

sub _truncate_name {
    my $name = shift;
    my @name_parts = split('\.', $name);
    return $name_parts[0];
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
