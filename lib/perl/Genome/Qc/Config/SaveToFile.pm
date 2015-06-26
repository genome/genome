package Genome::Qc::Config::SaveToFile;

use strict;
use warnings;
use Genome;
use JSON qw(decode_json);

class Genome::Qc::Config::SaveToFile {
    is => 'Command::V2',
    has => [
        config_item => {
            is => 'Genome::Qc::Config',
            doc => 'QC configuration to save to file',
        },
        file_path => {
            is => 'Path',
            doc => 'File path to store the QC configuration file',
        },
    ],
};

sub execute {
    my $self = shift;

    my $config = decode_json($self->config_item->config);
    Genome::Config::Parser::YAML->write($self->file_path, $config);

    return 1;
}

1;
