package Genome::Config::Translator;

use strict;
use warnings;

use Genome;

class Genome::Config::Translator {
    is => 'UR::Object',
    is_transactional => 0,
};

sub get_rule_model_map_from_config {
    my $class = shift;
    my $config = shift;

    my $config_hash = Genome::Config::Parser->parse($config->file_path);

    my @rules = Genome::Config::Rule->create_from_hash($config_hash->{rules});
    my $models = _wrap_single_model_hashes($config_hash->{models});
    $models = _add_config_profile_item_id($models, $config);

    return Genome::Config::RuleModelMap->create(
        rules => \@rules,
        models => $models,
        config => $config,
    );
}

sub _wrap_single_model_hashes {
    my $models = shift;

    while (my ($model_type, $model_config) = each %$models) {
        if (ref $model_config eq 'HASH') {
            $models->{$model_type} = [ $model_config ];
        }
    }

    return $models;
}

sub _add_config_profile_item_id {
    my $models = shift;
    my $config_profile_item = shift;

    while (my (undef, $model_config) = each %$models) {
        for (@{$model_config}) {
            $_->{config_profile_item} = $config_profile_item
        }
    }
    return $models;
}

1;
