package Genome::Config::AnalysisProject::Command::AddInstrumentDataToAnalysisProject;

use strict;
use warnings;

use Genome;

class Genome::Config::AnalysisProject::Command::AddInstrumentDataToAnalysisProject {
    is => 'Command::V2',
    has_input => [
       analysis_project  => {
            is                  => 'Genome::Config::AnalysisProject',
            doc                 => 'the analysis project to add the config file to',
            shell_args_position => 1,
        },
       instrument_data  => {
            is                  => 'Genome::InstrumentData::Imported',
            doc                 => 'imported instrument data to add to this analysis project',
            is_many             => 1,
            shell_args_position => 2,
        }
    ],
};

sub help_brief {
    return 'associate imported instrument data with an analysis project';
}

sub help_synopsis {
    return "genome config analysis-project add-instrument-data-to-analysis-project <analysis-project> <instrument-data>";
}

sub help_detail {
    return <<"EOS"
This will associate instrument data with an analysis project. It is limited to working on imported data only
as data produced in-house will have this information set up in LIMS
EOS
}

sub execute {
    my $self = shift;

    my @instrument_data = $self->instrument_data;
    my $analysis_project = $self->analysis_project;

    for my $instrument_data (@instrument_data) {
        Genome::Config::AnalysisProject::InstrumentDataBridge->create(
            analysis_project => $analysis_project,
            instrument_data => $instrument_data,
        );
    }

    return 1;
}

1;
