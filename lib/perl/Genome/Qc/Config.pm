package Genome::Qc::Config;

use strict;
use warnings;
use Genome;
use JSON qw(decode_json);
use Set::Scalar;

class Genome::Qc::Config {
    roles => [
        'Genome::Role::ObjectWithTimestamps',
        'Genome::Role::ObjectWithCreatedBy',
    ],
    table_name => 'config.qc',
    id_by => [
        id => {
            is => 'Text',
            len => 64,
        }
    ],
    has => [
        name => {
            is => 'Text',
        },
        type => {
            is => 'String',
            is_optional => 1,
            valid_values => ['wgs', 'exome', 'all'],
        },
        config => {
            is => 'Text',
        }
    ],
    schema_name => 'GMSchema',
    data_source => 'Genome::DataSource::GMSchema',
    id_generator => '-uuid',
};

sub __display_name__ { sprintf('%s for %s (%s)', $_[0]->name, $_[0]->type, $_[0]->id); }

sub __errors__ {
    my $self = shift;

    my @errors = $self->SUPER::__errors__(@_);

    my $available_tools = Set::Scalar->new(Genome::Qc::Tool->available_tools);
    my $config = decode_json($self->config);
    while ( my ($config_element, $tool_config) =  each %$config ) {
        for my $key ($self->required_keys_for_tool) {
            unless ($tool_config->{$key}) {
                push @errors, UR::Object::Tag->create(
                    type => 'error',
                    properties => ['config'],
                    desc => "Missing key '$key' for config element ($config_element)",
                );
            }
        }

        if (defined($tool_config->{class}) && !$available_tools->has($tool_config->{class})) {
            push @errors, UR::Object::Tag->create(
                type => 'error',
                properties => ['config'],
                desc => sprintf("Tool class (%s) for config element (%s) not found under Genome::Qc::Tool", $tool_config->{class}, $config_element),
            );
        }
    }

    return @errors;
}

sub required_keys_for_tool {
    return ('class', 'params');
}

sub get_commands_for_alignment_result {
    my $self = shift;
    my $is_capture = shift;

    my $config = decode_json($self->config);

    if ($is_capture && $self->_invalid_config_types_for_capture_data->contains($self->type)) {
        die $self->error_message('Config type (%s) of config (%s) is not valid for capture data', $self->type, $self->name);
    }
    if (!$is_capture && $self->_invalid_config_types_for_wgs_data->contains($self->type)) {
        die $self->error_message('Config type (%s) of config (%s) is not valid for whole genome data', $self->type, $self->name);
    }

    return $config;
}

sub _invalid_config_types_for_capture_data {
    return Set::Scalar->new('wgs');
}

sub _invalid_config_types_for_wgs_data {
    return Set::Scalar->new('exome');
}

1;

