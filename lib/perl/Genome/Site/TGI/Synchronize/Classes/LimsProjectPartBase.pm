package Genome::Site::TGI::Synchronize::Classes::LimsProjectPartBase; 

use strict;
use warnings;

use Genome;

class Genome::Site::TGI::Synchronize::Classes::LimsProjectPartBase {
    is => 'Genome::Site::TGI::Synchronize::Classes::LimsBase',
     schema_name => 'GMSchema',
    data_source => 'Genome::DataSource::Oltp',
};

sub entity_base_name {
    my $self = shift;
    my $class = $self->class;
    $class =~ s/Genome::Site::TGI::Synchronize::Classes::LimsProject//;
    return lc Genome::Utility::Text::camel_case_to_string($class);
}

sub entity_name {
    my $self = shift;
    return 'project '.$self->entity_base_name;
}

sub genome_class_for_comparison { 
    my $self = shift;
    return 'Genome::Site::TGI::Synchronize::Classes::'.join('', map { ucfirst } split(/\s+/, $self->entity_name));
}

sub genome_class_for_create { return 'Genome::ProjectPart'; }

sub properties_to_copy {
    my $self = shift;
    return ( 'project_id', 'entity_id' );
}

sub entity_class_name { 
    my $self = shift;
    return 'Genome::'.join('', map { ucfirst } split(/\s+/, $self->entity_base_name));
}

sub label{
    my $self = shift;
    join('_', split(/\s+/, $self->entity_base_name));
}

sub params_for_create_in_genome {
    my $self = shift;
    return (
        project_id => $self->project_id,
        entity_id => $self->entity_id,
        entity_class_name => $self->entity_class_name,
        label => $self->label,
    );
}

1;

