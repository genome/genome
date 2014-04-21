package Genome::Annotation::Factory;

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
    return ['Genome::Annotation', (map {'Genome::Annotation::' . $_} @{search_dirs()})];
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
    only => qr(Interpreter$),
    sub_name => 'interpreters';

use Module::Pluggable
    require => 1,
    search_path => search_path(),
    only => qr(Reporter$),
    sub_name => 'reporters';

class Genome::Annotation::Factory {};

sub names {
    my ($self, $accessor) = @_;
    return keys %{$self->_load($accessor)};
}

sub get_object {
    my ($self, $accessor, $name, $params) = validate_pos(@_, 1, 1, 1, 1);

    if (exists $self->_load($accessor)->{$name}) {
        my $pkg = $self->_load($accessor)->{$name};
        return $pkg->create(%{$params});
    } else {
        confess sprintf("No $accessor with name ($name) available $accessor are:\n    %s\n",
            join("\n    ", $self->names($accessor)));
    }
}

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

1;
