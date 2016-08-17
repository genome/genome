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
    my $anp = $self->analysis_project;


    for my $instrument_data ($anp->instrument_data) {
        my @allocations = Genome::Disk::Allocation->get(owner_id => $instrument_data->id, owner_class_name => $instrument_data->class);

        if (@allocations) {
            for my $allocation (@allocations) {
                if ($allocation->is_archived) {
                    say 'genome disk allocation unarchive --analysis-project ' . $anp->id . ' ' . $allocation->id;
                }
            }
        }

        $self->unarchive_additional_data($instrument_data);
    }
}

sub valid_statuses {
    return ('In Progress', 'Hold', 'Pending');
}

sub unarchive_additional_data {
    #this is a hook for overriding--by default do nothing;
    return 1;
}

1;
