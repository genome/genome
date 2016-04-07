package Genome::Qc::Command::Config::View;

use strict;
use warnings 'FATAL';

use Genome;

class Genome::Qc::Command::Config::View {
    is => 'Genome::Command::Viewer',
    has => {
        qc_config => {
            is => 'Genome::Qc::Config',
            shell_args_position => 1,
            doc => 'QC config to display.',
        },
    },
    doc => 'View a QC configuration',
};

sub help_detail { 'Displays basic information about a qc config' }

sub write_report {
    my ($self, $width, $handle) = @_;
    $self->_write_basic_info($width, $handle);
    $self->_write_decoded_config($width, $handle);
    return 1;
}

sub _write_basic_info {
    my ($self, $width, $handle) = @_;
    printf $handle "Name:   %s\n", $self->qc_config->name;
    printf $handle "ID:     %s\n", $self->qc_config->id;
    printf $handle "Type:   %s\n", $self->qc_config->type;
}

sub _write_decoded_config {
    my ($self, $width, $handle) = @_;
    my $string = $self->qc_config->config_to_yaml;
    $string =~ s/^-+/Config:/;
    print $handle $string;
}

1;

