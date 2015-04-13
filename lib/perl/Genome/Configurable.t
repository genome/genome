#!/usr/bin/env genome-perl

use strict;
use warnings;

use Genome qw();

use Test::More tests => 3;

use File::Temp qw();
use Genome::Test::Config qw(setup_config);
use Path::Class::Dir qw();

use_ok('Genome::Configurable');

my $temp_home_dir = File::Temp->newdir();
my $temp_conf_dir = File::Temp->newdir();
local $ENV{XGENOME_CONFIG_HOME} = $temp_home_dir->dirname;
local $ENV{XGENOME_CONFIG_DIRS} = $temp_conf_dir->dirname;

setup_config(
    home => {
        dir => Path::Class::Dir->new($temp_home_dir->dirname, 'genome'),
        config => {
            'foo.name' => 'bar',
        },
    },
    conf => [
        {
            dir => Path::Class::Dir->new($temp_conf_dir->dirname, 'genome'),
            spec => {
                'foo.name' => {
                    type => 'Str',
                },
            },
        },
    ],
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
