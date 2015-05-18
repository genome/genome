package Genome::Test::Config;

use strict;
use warnings;

use Exporter qw(import);
use File::Temp qw();
use Genome::Carp qw(croakf);
use Params::Validate qw(validate);
use Path::Class::Dir qw();
use Path::Class::File qw();
use POSIX qw(EEXIST);
use YAML::Syck qw();

our @EXPORT_OK = qw(setup_config);

sub setup_config {
    my %params = validate(@_, {
        spec   => { default => {} },
        home   => { default => {} },
        global => { default => {} },
    });

    for my $key (qw(XGENOME_CONFIG_SNAP XGENOME_CONFIG_HOME XGENOME_CONFIG_DIRS)) {
        unless ($ENV{$key}) {
            croakf '%s must be set before calling setup_config', $key;
        }
        unless (-d $ENV{$key}) {
            croakf '%s directory must be created before calling setup_config', $key;
        }
    }

    for my $key (keys %{$params{spec}}) {
        setup_spec_file(
            name => $key . '.yaml',
            data => $params{spec}->{$key},
        );
    }

    setup_home_file(data => $params{home});
    setup_global_file(data => $params{global});
}

sub setup_spec_file {
    my $dir = Path::Class::Dir->new($ENV{XGENOME_CONFIG_SNAP}, 'genome');
    write_yaml_data(dir => $dir, @_);
}

sub setup_home_file {
    my $dir = Path::Class::Dir->new($ENV{XGENOME_CONFIG_HOME}, 'genome');
    write_yaml_data(dir => $dir, name => 'config.yaml', @_);
}

sub setup_global_file {
    my $dir = Path::Class::Dir->new($ENV{XGENOME_CONFIG_DIRS}, 'genome');
    write_yaml_data(dir => $dir, name => 'config.yaml', @_);
}

sub write_yaml_data {
    my %params = validate(@_, {
        dir => 1,
        name => 1,
        data => 1,
    });

    $params{dir}->mkpath();
    my $config_file = Path::Class::File->new($params{dir}, $params{name});
    YAML::Syck::DumpFile($config_file, $params{data});
}

sub temp_dir_helper {
    my @temp_dirs;
    my $new_temp_dir = sub {
        push @temp_dirs, File::Temp->newdir();
        return $temp_dirs[-1]->dirname;
    };
    return (\@temp_dirs, $new_temp_dir);
}
