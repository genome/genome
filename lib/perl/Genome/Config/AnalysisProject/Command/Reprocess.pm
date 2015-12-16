package Genome::Config::AnalysisProject::Command::Reprocess;

use strict;
use warnings;

use Genome;

class Genome::Config::AnalysisProject::Command::Reprocess {
    is => 'Genome::Config::AnalysisProject::Command::Base',
    has_optional => {
        instrument_data => {
            is => 'Genome::InstrumentData',
            is_many => 1,
            doc => 'Only reprocess these instrument data.',
        },
    },
};

sub help_brief {
    return 'reschedules instrument data for processing';
}

sub help_synopsis {
    return 'genome config analysis-project reprocess <analysis-project>';
}

sub help_detail {
    return <<'EOS';
This will reschedule all instrument-data for the analysis-project to be
reprocessed.  This command has the same effect as the --reprocess-existing
option to the add-config-file and add-menu-item commands.

If the analysis-project is currently held or pending, the reprocessing will
not occur until the analysis-project has been released.
EOS
}

sub valid_statuses {
    return ('Pending', 'Hold', 'In Progress');
}

sub execute {
    my $self = shift;

    for my $bridge ( $self->_resolve_analysis_project_bridges ) {
        $bridge->reschedule;
    }

    return 1;
}

sub _resolve_analysis_project_bridges {
    my $self = shift;

    return $self->analysis_project->analysis_project_bridges if not $self->instrument_data;

    return Genome::Config::AnalysisProject::InstrumentDataBridge->get(
        analysis_project => $self->analysis_project,
        instrument_data => [ $self->instrument_data ],
    );
}

1
