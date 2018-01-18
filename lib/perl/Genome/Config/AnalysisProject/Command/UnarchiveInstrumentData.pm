package Genome::Config::AnalysisProject::Command::UnarchiveInstrumentData;

use strict;
use warnings;

use Genome;

use feature qw(say);

class Genome::Config::AnalysisProject::Command::UnarchiveInstrumentData {
    is => 'Genome::Config::AnalysisProject::Command::Base',
    has_input => [
        volume => {
            is => 'Genome::Disk::Volume',
            doc => 'A location to which scratch allocation information can be stored during unarchiving',
            is_optional => 1,
        },
    ],
    doc => 'a helper that produces commands that can be run to unarchive data',
};

sub execute {
    my $self = shift;

    my $unarchive_cmd = Genome::InstrumentData::Command::Unarchive->create(
        instrument_data => [$self->analysis_project->instrument_data],
        analysis_project => $self->analysis_project,
        volume => $self->volume,
    );
    unless ($unarchive_cmd) {
        $self->fatal_message('Failed to create unarchive instrument data command!');
    }
    unless ($unarchive_cmd->execute) {
        $self->fatal_message('Failed to execute unarchive instrument data command!');
    }

    return 1;
}

sub valid_statuses {
    return ('In Progress', 'Hold', 'Pending');
}

1;
