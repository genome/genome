package Genome::Site::TGI::Synchronize::Classes::LimsProjectPartBase; 

use strict;
use warnings;

use Genome;

class Genome::Site::TGI::Synchronize::Classes::LimsProjectPartBase {
    is => 'Genome::Site::TGI::Synchronize::Classes::LimsBase',
    has_constant_calculated => [
        entity_id => {
            calculate_from => [qw/ label /],
            calculate => q| my $method = $label.'_id'; return $self->$method; |, 
        },
        entity_class_name => { 
            calculate_from => [qw/ entity_base_name /],
            calculate => q| return 'Genome::'.join('', map { ucfirst } split(/\s+/, $entity_base_name)); |, 
        },
        label => {
            calculate_from => [qw/ entity_base_name /],
            calculate => q| join('_', split(/\s+/, $entity_base_name)); |, 
        },
    ],
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
    return (qw/ project_id entity_class_name entity_id label /);
}

1;

