package Genome::Qc::Config;

use strict;
use warnings;
use Genome;
use JSON qw(decode_json);
use Set::Scalar;

class Genome::Qc::Config {
    is => [
        'Genome::Utility::ObjectWithTimestamps',
        'Genome::Utility::ObjectWithCreatedBy',
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

sub get_commands_for_alignment_result {
    my $self = shift;
    my $is_capture = shift;

    my $config = decode_json($self->config);

    unless ($is_capture && $self->_config_types_for_capture_data->contains($self->type)) {
        die $self->error_message('Config type (%s) of config (%s) is not valid for capture data', $self->type, $self->name);
    }

    return $config;
}

sub _config_types_for_capture_data {
    return Set::Scalar->new('exome', 'all');
}

1;

