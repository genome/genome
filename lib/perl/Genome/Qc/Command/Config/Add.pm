package Genome::Qc::Command::Config::Add;

use strict;
use warnings;
use Genome;
use Set::Scalar;
use JSON qw(encode_json);

class Genome::Qc::Command::Config::Add {
    is => 'Command::V2',
    has => [
        name => {
            is => 'String',
            doc => 'Name of the qc configuration to add',
        },
        file_path => {
            is => 'Path',
            doc => 'File path to the qc configuration file to add',
        },
        type => {
            is => 'String',
            valid_values => ['wgs', 'exome', 'all'],
            doc => 'Type of data this qc configuration applies to',
        },
    ],
    doc => 'A command to add a QC configuration to the database',
};

sub help_detail {
    return 'A command to add a QC configuration to the database';
}

sub execute {
    my $self = shift;

    if (my @existing_configs = Genome::Qc::Config->get(name => $self->name)) {
        die $self->error_message("A config item with name (%s) already exists", $self->name);
    }

    my $config = Genome::Config::Parser::YAML->parse($self->file_path);
    my $config_item = Genome::Qc::Config->create(
        name => $self->name,
        type => $self->type,
        config => encode_json($config),
    );

    return 1;
}

1;
