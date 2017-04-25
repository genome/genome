package Genome::InstrumentData::Command::Unarchive;

use strict;
use warnings;

use Genome;

use feature qw(say);

class Genome::InstrumentData::Command::Unarchive {
    is => 'Command::V2',
    has_input => [
        instrument_data => {
            is => 'Genome::InstrumentData',
            doc => 'instrument data to unarchive',
            is_many => 1,
            shell_args_position => 1,
        },
        analysis_project => {
            is => 'Genome::Config::AnalysisProject',
            doc => 'the analysis project to assign the imported data to.',
            shell_args_position => 2,
        },
        volume => {
            is => 'Genome::Disk::Volume',
            doc => 'A location to which scratch allocation information can be stored during unarchiving',
            is_optional => 1,
        },
    ],
    doc => 'a helper that produces commands that can be run to unarchive instrument data and import the data into a GMS allocation',
};

sub execute {
    my $self = shift;

    for my $instrument_data ($self->instrument_data) {

        my @allocations = Genome::Disk::Allocation->get(owner_id => $instrument_data->id, owner_class_name => $instrument_data->class);

        if (@allocations) {
            for my $allocation (@allocations) {
                if ($allocation->is_archived) {
                    say 'genome disk allocation unarchive --analysis-project ' . $self->analysis_project->id . ' ' . $allocation->id;
                }
            }
        }
        $self->unarchive_additional_data($instrument_data);
    }
    return 1;
}

sub unarchive_additional_data {
    my $self = shift;
    
    #this is a hook for overriding--by default do nothing;
    $self->fatal_message('Failed to override unarchive_additional_data');
}

1;
