#!/usr/bin/env genome-perl

use strict;
use warnings;

use Genome qw();

use Test::More tests => 4;
use Test::Fatal qw(exception);

use File::Temp qw();
use Genome::Test::Config qw(setup_config);
use Path::Class::Dir qw();

use_ok('Genome::Config');

subtest 'basic lookup' => sub {
    plan tests => 2;

    my $temp_home_dir = File::Temp->newdir();
    my $temp_conf_dir = File::Temp->newdir();
    local $ENV{XGENOME_CONFIG_HOME} = $temp_home_dir->dirname;
    local $ENV{XGENOME_CONFIG_DIRS} = $temp_conf_dir->dirname;

    setup_config(
        home => {
            dir => Path::Class::Dir->new($temp_home_dir->dirname, 'genome'),
            config => {
                home_key => 'home_dir_value',
            },
        },
        conf => [
        {
            dir => Path::Class::Dir->new($temp_conf_dir->dirname, 'genome'),
            config => {
                conf_key => 'conf_dir_value',
                home_key => 'conf_dir_value',
            },
            spec => {
                home_key => {
                    type => 'Str',
                },
                conf_key => {
                    type => 'Str',
                },
            },
        },
        ],
    );

    is(Genome::Config::get('home_key'), 'home_dir_value', 'looked up correct value for home_key');
    is(Genome::Config::get('conf_key'), 'conf_dir_value', 'looked up correct value for conf_key');
};

subtest 'required value' => sub {
    plan tests => 2;

    my $temp_conf_dir = File::Temp->newdir();
    local $ENV{XGENOME_CONFIG_DIRS} = $temp_conf_dir->dirname;
    setup_config(
        home => {},
        conf => [
        {
            dir => Path::Class::Dir->new($temp_conf_dir->dirname, 'genome'),
            spec => {
                some_key => {
                    type => 'Str',
                },
                some_key_with_default => {
                    type => 'Str',
                    default_value => '',
                },
            },
        },
        ],
    );

    my $exception = exception { Genome::Config::get('some_key') };
    ok($exception, 'got exception');
    is(Genome::Config::get('some_key_with_default'), '', 'got empty string');
};

subtest 'validation' => sub {
    plan tests => 2;

    my $temp_conf_dir = File::Temp->newdir();
    local $ENV{XGENOME_CONFIG_DIRS} = $temp_conf_dir->dirname;
    setup_config(
        home => {},
        conf => [
        {
            dir => Path::Class::Dir->new($temp_conf_dir->dirname, 'genome'),
            config => {
                bad_numeric_key => 'abc',
                good_numeric_key => '123',
            },
            spec => {
                bad_numeric_key => {
                    type => 'Int',
                    validators => [ 'numeric' ],
                },
                good_numeric_key => {
                    type => 'Int',
                    validators => [ 'numeric' ],
                },
            },
        },
        ],
    );

    my $exception = exception { Genome::Config::get('bad_numeric_key') };
    ok($exception, 'got exception');
    is(Genome::Config::get('good_numeric_key'), '123', 'got number');
};
