package Genome::Timeline::Event::AnalysisProject;

use strict;
use warnings;

use Genome;

class Genome::Timeline::Event::AnalysisProject {
    is => 'Genome::Timeline::Event',
    table_name => 'timeline.analysis_project',
    id_by => [
        id => { is => 'Text', len => 64 },
    ],
    has => [
        object_class_name => {
            is => 'Text',
            is_constant => 1,
            default_value => 'Genome::Config::AnalysisProject',
            valid_values => ['Genome::Config::AnalysisProject'],
        },
        object_id => { is => 'Text' },
        analysis_project => {
            is => 'UR::Object',
            id_by => 'object_id',
            id_class_by => 'object_class_name',
            constraint_name => 'anp_event_anp_fk',
        },
        name => { is => 'Text' },
        type => {
            is => 'Genome::Timeline::Event::AnalysisProjectEventType',
            id_by => 'name',
            constraint_name => 'anp_event_type_fk',
        },
        status => {
            is => 'Text',
            len => 255,
            default_value => 'Pending',
            valid_values => [ "Pending", "Approved", "In Progress", "Completed", "Archived", "Hold" ],
        },
        is_cle => { is => 'Boolean', len => 1 },
        run_as => { is => 'Text', len => 64 },

        # added by "ur update classes-from-db" on 22 Sep 2014 - is this column used?
        reason => { is => 'Text' },
    ],
    schema_name => 'GMSchema',
    data_source => 'Genome::DataSource::GMSchema',
    id_generator => '-uuid',
};

Genome::Timeline::Event->_define_event_constructors(
    __PACKAGE__,
    Genome::Timeline::Event::AnalysisProjectEventType->__meta__->property('id')->valid_values(),
);

sub _properties_to_snapshot {
    return (
        'status' => 'status',
        'is_cle' => 'is_cle',
        'run_as' => 'run_as',
    );
}

1;
