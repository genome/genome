package Genome::Test::Config;

use strict;
use warnings;

use Exporter qw(import);
use Params::Validate qw(validate);
use Path::Class::File qw();
use YAML::Syck qw();

our @EXPORT_OK = qw(setup_config);

sub setup_config {
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
