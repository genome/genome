package Genome::Model::Tools::EdgeR::Base;

use Genome;
use Carp qw/confess/;

use strict;
use warnings;

class Genome::Model::Tools::EdgeR::Base {
    is => "Command::V2",
    has => [
        counts_file => {
            is => "Text",
            doc => "Tab delimited file of read counts per structure",
        },

        groups => {
            is => "Text",
            doc => "Comma separated list of conditions corresponding to each " .
                "column in counts_file (e.g, normal,tumor,...)",
        },

        output_file => {
            is => "Text",
            doc => "Output file path",
        },

        p_value => {
            is => "Number",
            doc => "The p-value used for significance testing (0 < p < 1)",
            default_value => 0.05,
        },
    ],
};

sub _validate_params {
    my $self = shift;

    if ($self->p_value <= 0 or $self->p_value >= 1) {
        confess "--p-value parameter must satisfy 0 < p < 1";
    }

    my @groups = split(",", $self->groups);
    my %group_sizes;
    for my $g (@groups) {
        return if ++$group_sizes{$g} > 1;
    }

    confess sprintf("At least one group must contain replication! Groups " .
        "assignments were: %s\n",
        $self->groups);
}

sub execute {
    my $self = shift;

    $self->_validate_params;

    my $cmd = $self->construct_r_command;

    return Genome::Sys->shellcmd(
        cmd => $cmd
        );
}

sub construct_r_command {
    my $self = shift;

    $self->error_message("R command must be specified in the child class");

    return 0;
}

1;
