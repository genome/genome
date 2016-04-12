package Genome::Qc::Command::Config::SaveToFile;

use strict;
use warnings;
use Genome;
use JSON qw(decode_json);

class Genome::Qc::Command::Config::SaveToFile {
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
    doc => 'A command to write an existing QC configuration to a file',
};

sub help_detail {
    return 'A command to write an existing QC configuration to a file';
}

sub execute {
    my $self = shift;
    return Genome::Sys->write_file($self->file_path, $self->config_item->config_to_yaml);
}

1;
