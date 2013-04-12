package Genome::Config::MaskedConfigurationReader;
use warnings;
use strict;

use Genome;

class Genome::Config::MaskedConfigurationReader {
    is => 'UR::Object',
    has => [
        config_handler => {
            is => 'UR::Object',
        },
        mask_handler => {
            is => 'UR::Object',
        },
        default_handler => {
            is => 'UR::Object',
        },
        configuration_parser => {
            is => 'UR::Object',
        },
        configuration_copy_strategy => {
            is => 'UR::Object',
        },
    ],
    doc => 'takes a configuration mask and a configuration set and allows it to be queried',
};

sub get_config {
    my ($self, %params) = @_;
    my @params = values %params;
    return unless $self->mask_handler->parameters_exist(@params);

    my $existing_config = $self->config_handler->parameters_exist(@params);
    if ($existing_config && -f $existing_config){
        return $self->configuration_parser->parse($existing_config);
    } else {
        my $config_value = $self->default_handler->get_leaf(@params);
        $self->configuration_copy_strategy->copy_config($config_value,
            $self->default_handler->base_path, $self->config_handler->base_path);
        return $self->configuration_parser->parse($config_value);
    }
}

1;
