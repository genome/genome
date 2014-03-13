package Genome::Timeline::Event::AnalysisProject;

use strict;
use warnings;

use Genome;

class Genome::Timeline::Event::AnalysisProject {
    is => 'Genome::Timeline::Event',
    table_name => 'timeline.analysis_project',
    has => [
        object_class_name => {
            is_constant => 1,
            default_value => 'Genome::Config::AnalysisProject',
            valid_values => ['Genome::Config::AnalysisProject'],
        },
        analysis_project => {
            id_by => 'object_id',
            is => 'UR::Object',
            id_class_by => 'object_class_name'
        },
        type => {
            is => 'Genome::Timeline::Event::AnalysisProjectEventType',
            id_by => 'name',
        },
        status => {
            is => 'Text',
            default_value => 'Pending',
            valid_values => ['Pending', 'Approved', 'In Progress', 'Completed', 'Archived', 'Hold'],
        },
        is_cle => {
            is => 'Boolean',
        },
        run_as => {
            is => 'Text',
        },
    ],
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
