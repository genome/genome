#!/usr/bin/env genome-perl

use strict;
use warnings;

use Test::More tests => 2;

use File::Temp qw();
use Params::Validate qw(validate);
use Path::Class qw();
use YAML::Syck qw();

use_ok('Genome::Config');

subtest 'get: basic' => sub {
    plan tests => 2;

    my $temp_home_dir = File::Temp->newdir();
    my $temp_conf_dir = File::Temp->newdir();
    local $ENV{XGENOME_ENABLE_USER_CONFIG} = 1;
    local $ENV{XGENOME_CONFIG_HOME} = $temp_home_dir->dirname;
    local $ENV{XGENOME_CONFIG_DIRS} = $temp_conf_dir->dirname;

    setup(
        home => {
            dir => Path::Class::Dir->new($temp_home_dir->dirname, 'genome'),
            config => {
                home_key => 'home_dir_value',
            },
            spec => {
                home_key => {
                    type => 'Str',
                },
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

sub setup {
    my %params = validate(@_, {
        home => 1,
        conf => 1,
    });

    for my $blob ($params{home}, @{$params{conf}}) {
        unless (-d $blob->{dir}) {
            mkdir $blob->{dir};
        }
        for my $key (keys %{$blob->{spec}}) {
            setup_spec_file(
                dir => $blob->{dir},
                key => $key,
                spec => $blob->{spec}->{$key}
            );
        }
        setup_config_file(
            dir => $blob->{dir},
            data => $blob->{config},
        );
    }
}

sub setup_spec_file {
    my %params = validate(@_, {
        dir => 1,
        key => 1,
        spec => 1,
    });
    my $spec_file = Path::Class::File->new($params{dir}, $params{key} . '.yaml');
    YAML::Syck::DumpFile($spec_file . '', $params{spec});
}

sub setup_config_file {
    my %params = validate(@_, {
        dir => 1,
        data => 1,
    });
    my $config_file = Path::Class::File->new($params{dir}, 'config.yaml');
    YAML::Syck::DumpFile($config_file . '', $params{data});
}
