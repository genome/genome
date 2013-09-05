package Genome::Config::AnalysisProject::SubjectPairing;

use strict;
use warnings;

use Genome;

class Genome::Config::AnalysisProject::SubjectPairing {
    is           => ['Genome::Utility::ObjectWithTimestamps', 'Genome::Utility::ObjectWithCreatedBy'],
    id_generator => '-uuid',
    data_source  => 'Genome::DataSource::GMSchema',
    table_name   => 'subject.pairing',
    id_by        => [
        id => {
            is => 'Text',
            len => 64,
        },
    ],
    has => [
        analysis_project => {
            is    => 'Genome::Config::AnalysisProject',
            id_by => 'analysis_project_id'
        },
        control_subject => {
            is    => 'Genome::Subject',
            id_by => 'control_subject_id',
        },
        experimental_subject => {
            is    => 'Genome::Subject',
            id_by => 'experimental_subject_id',
        },
    ]
};

sub __display_name__ {
    my $self = shift;
    return sprintf('Control: %s, Experimental: %s as part of %s',
        $self->control_subject->__display_name__,
        $self->experimental_subject->__display_name__,
        $self->analysis_project->__display_name__);
}


1;
