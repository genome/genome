#!/usr/bin/env genome-perl

use strict;
use warnings;

use Genome qw();
use Genome::Test::Config qw(setup_config);

use Test::More tests => 3;

use_ok('Genome::Configurable');

my ($temp_dirs, $new_temp_dir) = Genome::Test::Config::temp_dir_helper();
local $ENV{XGENOME_CONFIG_SNAP} = $new_temp_dir->();
local $ENV{XGENOME_CONFIG_HOME} = $new_temp_dir->();
local $ENV{XGENOME_CONFIG_DIRS} = $new_temp_dir->();
setup_config(
    spec => {
        'foo.name' => {
            type => 'Str',
        },
    },
    home => {
        'foo.name' => 'bar',
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

my $foo1 = Genome::Foo->create();
is($foo1->name, 'bar');

my $foo2 = Genome::Foo->create(name => 'Joe');
is($foo2->name, 'Joe');
