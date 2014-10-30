package Genome::Config::AnalysisProject::Command::AddInstrumentData;

use strict;
use warnings;

use Genome;
use Set::Scalar;

class Genome::Config::AnalysisProject::Command::AddInstrumentData {
    is => 'Command::V2',
    has_input => [
       analysis_project => {
            is                  => 'Genome::Config::AnalysisProject',
            doc                 => 'the analysis project to which to add the instrument data',
            shell_args_position => 1,
       },
       instrument_data  => {
            is                  => 'Genome::InstrumentData::Imported',
            doc                 => 'imported instrument data to add to this analysis project',
            is_many             => 1,
            shell_args_position => 2,
       },
    ],
};

sub help_brief {
    return 'associate imported instrument data with an analysis project';
}

sub help_synopsis {
    return "genome config analysis-project add-instrument-data <analysis-project> <instrument-data>";
}

sub help_detail {
    return <<"EOS"
This will associate instrument data with an analysis project. It is limited to working on imported data only
as data produced in-house will have this information set up in LIMS
EOS
}

sub execute {
    my $self = shift;

    my $analysis_project = $self->analysis_project;
    my $already_assigned = Set::Scalar->new($analysis_project->instrument_data);
    my $desired_assigned = Set::Scalar->new($self->instrument_data);

    my $should_assign = $desired_assigned - $already_assigned;
    my $previously_assigned = $already_assigned->intersection($desired_assigned);

    $self->status_message('Asked to asssign (%s) of which (%s) were already assigned. Proceeding to assign (%s)',
        $desired_assigned->size, $previously_assigned->size, $should_assign->size);

    for my $instrument_data ($should_assign->members) {
        Genome::Config::AnalysisProject::InstrumentDataBridge->create(
            analysis_project => $analysis_project,
            instrument_data => $instrument_data,
        );
    }

    return 1;
}

1;
