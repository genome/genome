package Genome::Config::AnalysisProject::SubjectMapping;

use strict;
use warnings;

use Genome;

class Genome::Config::AnalysisProject::SubjectMapping {
    roles => ['Genome::Role::ObjectWithTimestamps', 'Genome::Role::ObjectWithCreatedBy'],
    table_name => 'config.subject_mapping',
    id_by => [
        id => { is => 'Text', len => 64 },
    ],
    has => [
        analysis_project => {
            is => 'Genome::Config::AnalysisProject',
            id_by => 'analysis_project_id',
            constraint_name => 'subject_mapping_analysis_project_id_fkey',
        },
        subject_bridges => {
            is => 'Genome::Config::AnalysisProject::SubjectMapping::Subject',
            reverse_as => 'subject_mapping',
            is_many => 1,
        },
        subjects => {
            via => 'subject_bridges',
            to => 'subject',
        },
        inputs => {
            is => 'Genome::Config::AnalysisProject::SubjectMapping::Input',
            reverse_as => 'subject_mapping',
            is_many => 1,
        },
        tag_bridges => {
            is => 'Genome::Config::Tag::AnalysisProject::SubjectMapping',
            reverse_as => 'subject_mapping',
            is_many => 1,
        },
        tags => {
            is => 'Genome::Config::Tag',
            via => 'tag_bridges',
            to => 'tag',
            is_many => 1,
            is_mutable => 1,
        },
    ],
    schema_name => 'GMSchema',
    data_source => 'Genome::DataSource::GMSchema',
    id_generator => '-uuid',
};

sub delete {
    my $self = shift;
    eval {
        for ($self->subject_bridges) {
            $_->delete();
        }
        for ($self->inputs) {
            $_->delete();
        }
        for ($self->tag_bridges) {
            $_->delete();
        }
    };
    if(my $error = $@) {
        die($error);
    }
    return $self->SUPER::delete();
}

1;
