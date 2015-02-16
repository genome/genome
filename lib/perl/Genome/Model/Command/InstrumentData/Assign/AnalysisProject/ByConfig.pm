package Genome::Model::Command::InstrumentData::Assign::AnalysisProject::ByConfig;

use strict;
use warnings;

use Genome;

class Genome::Model::Command::InstrumentData::Assign::AnalysisProject::ByConfig {
    is => 'Genome::Model::Command::InstrumentData::Assign::AnalysisProject',
    has_input => [
        config_profile_item => {
            is => 'Genome::Config::Profile::Item',
        },
    ],
};

sub _resolve_instrument_data {
    my $self = shift;
    my $model = $self->model;
    my $ruleset = Genome::Config::Translator->get_rule_model_map_from_config($self->config_profile_item);
    return grep { $ruleset->match($_) } $self->SUPER::_resolve_instrument_data();
}

1;
