#!/usr/bin/env genome-perl

use strict;
use warnings;

use Genome qw();

use Test::More tests => 2;

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
