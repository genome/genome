package Genome::Annotation::Factory;

use strict;
use warnings FATAL => 'all';
use Genome;
use Carp qw(confess);
use Memoize qw();
use File::Spec;
use File::Basename qw(dirname);
use Cwd qw(abs_path);

sub search_path {
    return [map {'Genome::Annotation::' . $_} @{search_dirs()}];
}

sub search_dirs {
    my @dirs;
    my $this_dir = dirname(abs_path(__FILE__));
    for my $candidate (glob(File::Spec->join($this_dir,'*'))) {
        push @dirs, $candidate if -d $candidate;
    }
    return [map {File::Spec->abs2rel($_, $this_dir)} @dirs];
}

use Module::Pluggable
    require => 1,
    search_path => search_path(),
    only => qr(Expert$),
    sub_name => 'experts';

use Module::Pluggable
    require => 1,
    search_path => search_path(),
    only => qr(Filter$),
    sub_name => 'filters';

use Module::Pluggable
    require => 1,
    search_path => search_path(),
    only => qr(Gatherer$),
    sub_name => 'gatherers';

use Module::Pluggable
    require => 1,
    search_path => search_path(),
    only => qr(View$),
    sub_name => 'views';

class Genome::Annotation::Factory {};

sub expert_names {
    my $self = shift;
    return $self->_names('experts');
}

sub get_expert {
    my $self = shift;
    my $name = shift;
    my $class = $self->_get('experts', $name);
    return $class->create(@_);
}

sub filter_names {
    my $self = shift;
    return $self->_names('filters');
}

sub get_filter {
    my $self = shift;
    my $name = shift;
    my $class = $self->_get('filters', $name);
    return $class->create(@_);
}

sub gatherer_names {
    my $self = shift;
    return $self->_names('gatherers');
}

sub get_gatherer {
    my $self = shift;
    my $name = shift;
    my $class = $self->_get('gatherers', $name);
    return $class->create(@_);
}

sub view_names {
    my $self = shift;
    return $self->_names('views');
}

sub get_view {
    my $self = shift;
    my $name = shift;
    my $class = $self->_get('views', $name);
    return $class->create(@_);
}

sub _names {
    my ($self, $accessor) = @_;
    return keys %{$self->_load($accessor)};
}

sub _get {
    my ($self, $accessor, $name) = @_;

    if (exists $self->_load($accessor)->{$name}) {
        return $self->_load($accessor)->{$name};
    } else {
        confess sprintf("No $accessor with name ($name) available $accessor are:\n    %s",
            join("\n    ", $self->_names($accessor)));
    }
}
Memoize::memoize('_get');

sub _load {
    my ($self, $accessor) = @_;

    my %plugins;
    for my $plugin ($self->$accessor) {
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
Memoize::memoize('_load');

1;
