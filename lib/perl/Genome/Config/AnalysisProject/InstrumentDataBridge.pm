package Genome::Config::AnalysisProject::InstrumentDataBridge;

use strict;
use warnings;

use Genome;

class Genome::Config::AnalysisProject::InstrumentDataBridge {
    is => 'Genome::Utility::ObjectWithTimestamps',
    table_name => 'config.instrument_data_analysis_project_bridge',
    id_by => [
        id => { is => 'Text', len => 64 },
    ],
    has => [
        analysis_project => {
            is => 'Genome::Config::AnalysisProject',
            id_by => 'analysis_project_id',
            constraint_name => 'instrument_data_analysis_project_bridg_analysis_project_id_fkey',
        },
        instrument_data => {
            is => 'Genome::InstrumentData',
            id_by => 'instrument_data_id',
            constraint_name => 'instrument_data_analysis_project_bridge_instrument_data_id_fkey',
        },
        status => {
            is => 'Text',
            default_value => 'new',
            valid_values => [ "new", "failed", "processed", "skipped" ],
        },
        fail_count => { is => 'Integer', len => 4, default_value => 0 },
        reason => { is => 'Text', is_optional => 1 },
    ],
    schema_name => 'GMSchema',
    data_source => 'Genome::DataSource::GMSchema',
    id_generator => '-uuid',
};

sub reschedule {
    my $self = shift;

    $self->status('new');
}

1;
