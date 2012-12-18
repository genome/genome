package Genome::Model::Build::MetagenomicComposition16s::GenerateReports;

use strict;
use warnings;

use Genome;

class Genome::Model::Build::MetagenomicComposition16s::GenerateReports {
    is => 'Command::V2',
    has_input => [
        build => { 
            is => 'Genome::Model::Build::MetagenomicComposition16s',
            is_output => 1,
        },
    ],
};

sub execute {
    my $self = shift;

    unless ( $self->build->generate_reports ) {
        $self->error_message("Failed to generate reports for ".$self->build->description);
        return;
    }

    return 1;
}

1;

