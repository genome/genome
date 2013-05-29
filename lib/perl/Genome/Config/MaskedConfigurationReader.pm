package Genome::Config::MaskedConfigurationReader;
use warnings;
use strict;

use Genome;

class Genome::Config::MaskedConfigurationReader {
    is => 'UR::Object',
    has => [
        config_handler => {
            is => 'Genome::Config::Handler',
        },
        mask_handler => {
            is => 'Genome::Config::Handler',
        },
        default_handler => {
            is => 'Genome::Config::Handler',
        },
        configuration_parser => {
            is => 'Genome::Config::Parser',
        },
        configuration_copy_strategy => {
            is => 'Genome::Config::CopyStrategy',
        },
    ],
    doc => 'takes a configuration mask and a configuration set and allows it to be queried',
};

sub get_config {
    my ($self, %params) = @_;
    my @params = values %params;
    return unless $self->mask_handler->valid_params(@params);

    my $existing_config = $self->config_handler->valid_params(@params);
    if ($existing_config && -f $existing_config){
        return $self->configuration_parser->parse($existing_config);
    } else {
        my $config_value = $self->default_handler->get_config(@params);
        $self->configuration_copy_strategy->copy_config($config_value,
            $self->default_handler->base_path, $self->config_handler->base_path);
        return $self->configuration_parser->parse($config_value);
    }
}

1;
