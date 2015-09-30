package Genome::Model::Tools::Analysis::RunWgsQc;

use strict;
use warnings;

use Genome;
use Command;
use IO::File;
use Sort::Naturally qw( nsort );

class Genome::Model::Tools::Analysis::RunWgsQc {
    is => 'Command',
    has => [
    model_group => {
        type => 'String',
        is_optional => 1,
        doc => "The ID of the Alignment Model Group to use",
    },
    ]
};

sub help_brief {
    "Extracts quality control metrics desired by production for a set of WGS (X Ten) models"
}


sub help_detail {
    <<'HELP';
HELP
}


sub execute {
    my $self=shift;

    warn "Extracting QC metrics for " . $self->model_group . "\n";


    ## The body of the program goes here ##
    




    ## Don't remove the return line ##
    return 1;
}


1;


