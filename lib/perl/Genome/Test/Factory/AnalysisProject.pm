package Genome::Test::Factory::AnalysisProject;
use base qw(Genome::Test::Factory::Base);

use strict;
use warnings;
use Genome;

use Genome::Test::Factory::Config::Profile::Item;

our @required_params = qw(name);

sub generate_obj {
    my $self = shift;
    my %params = @_;
    my $config_hash = delete $params{'config_hash'};

    unless(exists $params{status}) {
        $params{status} = 'In Progress'; #most commonly needed for tests
    }
    unless(exists $params{run_as}) {
        $params{run_as} = Genome::Sys->username;
    }

    my $project = Genome::Config::AnalysisProject->create(%params);

    if ($config_hash) {
        my @configs = (ref $config_hash eq 'ARRAY'? @$config_hash : $config_hash);
        for my $config (@configs) {
            my $config_profile_item = Genome::Test::Factory::Config::Profile::Item->setup_object(analysis_project => $project);
            $config->{config_profile_item} = $config_profile_item;
            $config_profile_item->{__dummy_config_hash__} = $config;
        }
    } else {
        Genome::Test::Factory::Config::Profile::Item->setup_object(analysis_project => $project);
    }

    return $project;

}

sub create_name {
    return Genome::Test::Factory::Util::generate_name('analysis_project_name');
}

use Genome::Config::Profile::Item;
use Genome::Config::Parser;
my $old_file_path_method = \&Genome::Config::Profile::Item::file_path;
my $old_is_concrete_method = \&Genome::Config::Profile::Item::is_concrete;
my $old_parse_method = \&Genome::Config::Parser::parse;
no warnings qw(redefine);

my $next = 0;
my %configs;

*Genome::Config::Profile::Item::file_path = sub {
    my $self = shift;
    if ($self->{__dummy_config_hash__}) {
        unless ($self->{__dummy_key__}) {
            my $key = '!dummy' . $next++;
            $configs{$key} = $self->{__dummy_config_hash__};
            $self->{__dummy_key__} = $key;
        }
        return $self->{__dummy_key__};
    } else {
        return $old_file_path_method->($self, @_);
    }
};
*Genome::Config::Profile::Item::is_concrete = sub {
    my $self = shift;
    return 1 if $self->{__dummy_config_hash__};
    return $old_is_concrete_method->($self, @_);
};
*Genome::Config::Parser::parse = sub {
    my $class = shift;
    my $path = shift;
    if (substr($path,0,1) eq '!') {
        return $configs{$path};
    } else {
        return $old_parse_method($class, $path, @_);
    }
};
use warnings;

1;
