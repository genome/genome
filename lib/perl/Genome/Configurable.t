#!/usr/bin/env genome-perl

use strict;
use warnings;

use Genome qw();
use Genome::Test::Config qw(setup_config);

use Test::More tests => 5;

use_ok('Genome::Configurable');

my ($temp_dirs, $new_temp_dir) = Genome::Test::Config::temp_dir_helper();
local $ENV{XGENOME_CONFIG_SNAP} = $new_temp_dir->();
local $ENV{XGENOME_CONFIG_HOME} = $new_temp_dir->();
local $ENV{XGENOME_CONFIG_DIRS} = $new_temp_dir->();
setup_config(
    spec => {
        'foo.name' => {},
        'bar.name' => {
            env => 'XGENOME_BAR',
        },
    },
    home => {
        'foo.name' => 'foo_value',
    },
);

UR::Object::Type->define(
    class_name => 'Genome::Foo',
    is => ['Genome::Configurable'],
    has => [
        name => {
            config => 'foo.name',
        },
    ],
);

UR::Object::Type->define(
    class_name => 'Genome::Bar',
    is => ['Genome::Configurable'],
    has => [
        name => {
            config => [qw(bar.name foo.name)],
        },
    ],
);

my $foo1 = Genome::Foo->create();
is($foo1->name, 'foo_value');

my $foo2 = Genome::Foo->create(name => 'Joe');
is($foo2->name, 'Joe');

{
    my $bar = Genome::Bar->create();
    is($bar->name, 'foo_value', q(got 'foo_value' when 'bar.name' is not set));
}

{
    my $guard = Genome::Config::set_env('bar.name', 'bar_value');
    my $bar = Genome::Bar->create();
    is($bar->name, 'bar_value', q(got 'bar_value' when 'bar.name' was set));
}
