package Genome::Qc::Config::Add;

use strict;
use warnings;
use Genome;
use Set::Scalar;
use JSON qw(encode_json);

class Genome::Qc::Config::Add {
    is => 'Command::V2',
    has => [
        name => {
            is => 'String',
        },
        file_path => {
            is => 'Path',
        },
        type => {
            is => 'String',
            valid_values => ['wgs', 'exome', 'all'],
        },
    ],
};

sub execute {
    my $self = shift;

    if (my @existing_configs = Genome::Qc::Config->get(name => $self->name)) {
        die $self->error_message("A config item with name (%s) already exists", $self->name);
    }

    my $config = Genome::Config::Parser::YAML->parse($self->file_path);
    $self->_validate_config($config);

    my $config_item = Genome::Qc::Config->create(
        name => $self->name,
        type => $self->type,
        config => encode_json($config),
    );

    return 1;
}

sub _validate_config {
    my $self = shift;
    my $config = shift;

    my $available_tools = Set::Scalar->new(Genome::Qc::Tool->available_tools);

    while ( my ($config_element, $tool_config) =  each %$config ) {
        for my $key ($self->required_keys_for_tool) {
            unless ($tool_config->{$key}) {
                die $self->error_message ("Missing key '$key' for config element (%s)", $config_element);
            }
        }

        unless ($available_tools->has($tool_config->{class})) {
            die $self->error_message ("Tool class (%s) for config element (%s) not found under Genome::Qc::Tool", $tool_config->{class}, $config_element);
        }
    }
}

sub required_keys_for_tool {
    return ('class', 'params');
}

1;
