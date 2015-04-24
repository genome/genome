#!/usr/bin/env genome-perl

use strict;
use warnings;

use Genome qw();
use Genome::Test::Config qw(setup_config);

use Test::More tests => 6;
use Test::Fatal qw(exception);

use_ok('Genome::Config');

subtest 'basic lookup' => sub {
    plan tests => 2;

    my ($temp_dirs, $new_temp_dir) = Genome::Test::Config::temp_dir_helper();
    local $ENV{XGENOME_CONFIG_SNAP} = $new_temp_dir->();
    local $ENV{XGENOME_CONFIG_HOME} = $new_temp_dir->();
    local $ENV{XGENOME_CONFIG_DIRS} = $new_temp_dir->();
    setup_config(
        spec => {
            home_key => {},
            conf_key => {},
        },
        home => {
            home_key => 'home_dir_value',
        },
        global => {
            conf_key => 'conf_dir_value',
            home_key => 'conf_dir_value',
        },
    );

    is(Genome::Config::get('home_key'), 'home_dir_value', 'looked up correct value for home_key');
    is(Genome::Config::get('conf_key'), 'conf_dir_value', 'looked up correct value for conf_key');
};

subtest 'required value' => sub {
    plan tests => 2;

    my ($temp_dirs, $new_temp_dir) = Genome::Test::Config::temp_dir_helper();
    local $ENV{XGENOME_CONFIG_SNAP} = $new_temp_dir->();
    local $ENV{XGENOME_CONFIG_HOME} = $new_temp_dir->();
    local $ENV{XGENOME_CONFIG_DIRS} = $new_temp_dir->();
    setup_config(
        spec => {
            some_key => {},
            some_key_with_default => {
                default_value => '',
            },
        },
    );

    my $exception = exception { Genome::Config::get('some_key') };
    ok($exception, 'got exception');
    is(Genome::Config::get('some_key_with_default'), '', 'got empty string');
};

subtest 'validation' => sub {
    plan tests => 2;

    my ($temp_dirs, $new_temp_dir) = Genome::Test::Config::temp_dir_helper();
    local $ENV{XGENOME_CONFIG_SNAP} = $new_temp_dir->();
    local $ENV{XGENOME_CONFIG_HOME} = $new_temp_dir->();
    local $ENV{XGENOME_CONFIG_DIRS} = $new_temp_dir->();
    setup_config(
        spec => {
            bad_numeric_key => {
                validators => [ 'numeric' ],
            },
            good_numeric_key => {
                validators => [ 'numeric' ],
            },
        },
        global => {
            bad_numeric_key => 'abc',
            good_numeric_key => '123',
        },
    );

    my $exception = exception { Genome::Config::get('bad_numeric_key') };
    ok($exception, 'got exception');
    is(Genome::Config::get('good_numeric_key'), '123', 'got number');
};

subtest 'sticky' => sub {
    plan tests => 2;

    my ($temp_dirs, $new_temp_dir) = Genome::Test::Config::temp_dir_helper();
    local $ENV{XGENOME_CONFIG_SNAP} = $new_temp_dir->();
    local $ENV{XGENOME_CONFIG_HOME} = $new_temp_dir->();
    local $ENV{XGENOME_CONFIG_DIRS} = $new_temp_dir->();
    setup_config(
        spec => {
            bad_sticky_key => {
                sticky => 1,
            },
            good_sticky_key => {
                sticky => 1,
                env => 'GOOD_STICKY_KEY',
            },
        },
        global => {
            bad_sticky_key => 'abc',
            good_sticky_key => 'abc',
        },
    );

    my $exception = exception { Genome::Config::get('bad_sticky_key') };
    ok($exception, 'got exception');
    is(Genome::Config::get('good_sticky_key'), 'abc', 'got value');
};

subtest 'env' => sub {
    plan tests => 2;

    my ($temp_dirs, $new_temp_dir) = Genome::Test::Config::temp_dir_helper();
    local $ENV{XGENOME_CONFIG_SNAP} = $new_temp_dir->();
    local $ENV{XGENOME_CONFIG_HOME} = $new_temp_dir->();
    local $ENV{XGENOME_CONFIG_DIRS} = $new_temp_dir->();
    setup_config(
        spec => {
            some_key => {
                env => 'SOME_KEY',
            },
        },
    );

    my $exception = exception { Genome::Config::get('some_key') };
    ok($exception, 'got exception');
    local $ENV{SOME_KEY} = 'abc';
    is(Genome::Config::get('some_key'), 'abc', 'got value');
};
