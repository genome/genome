package Genome::Config::Profile;

use strict;
use warnings;

use Genome;
use UR::Util;

use feature qw(switch);

class Genome::Config::Profile {
    is => 'UR::Object',
    is_transactional => 0,
    has => [
        config_rule_maps => {
            is_many => 1,
            is => 'Genome::Config::RuleModelMap'
        },
        analysis_project => {
            is => 'Genome::Config::AnalysisProject',
        },
    ]
};


sub create_from_analysis_project {
    my $class = shift;
    my $analysis_project = shift;

    my @config_rule_maps = map {
        Genome::Config::Translator->get_rule_model_map_from_config($_)
    } $analysis_project->config_items;

    return $class->create(
        config_rule_maps => \@config_rule_maps,
        analysis_project => $analysis_project,
    );
}

sub get_config {
    my $self = shift;
    my $inst_data = shift;

    my @model_hashes = map { $_->models }
        grep { $_->match($inst_data) }
        $self->config_rule_maps;

    return $self->_merge_extra_parameters($self->_merge_model_hashes(@model_hashes));
}

sub _merge_model_hashes {
    my $self = shift;
    my @source_hashes = @_;

    my $destination_hash = {};

    for my $source (@source_hashes) {
        for my $key (keys %$source) {
            $destination_hash->{$key} ||= [];
            push @{$destination_hash->{$key}}, @{$source->{$key}};
        }
    }

    return UR::Util::deep_copy($destination_hash);
}

sub _merge_extra_parameters {
    my $self = shift;
    my $config_hash = shift;

    for my $model_config (values %$config_hash) {
        for (@$model_config) {
             $self->_add_user_if_not_present($_);
        }
    }

    return $config_hash;
}

sub _add_user_if_not_present {
    my $self = shift;
    my $hashref = shift;

    $hashref->{user_name} ||= $self->analysis_project->run_as;

    return 1;
}

1;
