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

1;