package Genome::Site::TGI::Observers::ModelGroup;

use strict;
use warnings;

our %projects_with_deleted_model_groups;
our %deleted_bridges;

Genome::ModelGroup->add_observer(
    aspect => 'create',
    callback => \&model_group_create,
);

Genome::ModelGroup->add_observer(
    aspect => 'delete',
    callback => \&model_group_delete,
);

Genome::ModelGroup->add_observer(
    aspect => 'name',
    callback => \&model_group_rename,
);

Genome::ModelGroupBridge->add_observer(
    aspect => 'create',
    callback => \&model_group_bridge_create,
);

Genome::ModelGroupBridge->add_observer(
    aspect => 'delete',
    callback => \&model_group_bridge_delete,
);

sub model_group_create {
    my $self = shift;

    my $project = $self->project;
    return 1 if $project;

    $self->status_message('Create associated project');
    my $uuid = Genome::Project->__meta__->autogenerate_new_object_id_uuid;
    if ( not $uuid ) {
        die 'Failed to get uuid from Genome::Project! Cannot create model group.';
    }
    $self->uuid($uuid);
    my $name = $self->name;
    $project = Genome::Project->create(
        id => $uuid,
        name => $name,
    );
    if ( not $project ) {
        die 'Failed to create project to match model group.';
    }
    $self->name( $project->name ) if $project->name ne $self->name;
    $self->user_name( $project->parts(role => 'creator')->entity->email );

    for my $model ( $self->models ) {
        $project->add_part(entity => $model);
    }

    return 1;
}

sub model_group_delete {
    my $self = shift;

    my $project_id = $self->uuid;
    return 1 unless $project_id;
    $projects_with_deleted_model_groups{$project_id}++;

    return 1 if $Genome::Site::TGI::Observers::Project::groups_with_deleted_projects{$project_id};

    my $project = $self->project;
    return 1 if not $project;

    $self->status_message('Deleting associated project: '.$project->id);

    if ( not $project->delete ) {
        die 'Failed to delete project: '.$project->id;
    }

    return 1;
}

sub model_group_rename {
    my ($self, $property_name, $old_name, $new_name) = @_;

    if ( $self->uuid and my $project = Genome::Project->get($self->uuid) ) {
        if ( $project and $project->name ne $new_name ) {
            $self->status_message('Rename associated project');
            $project->rename($new_name);
        }
    }
    
    return 1 if $self->name eq $new_name;

    $self->rename($new_name) or die "Failed to rename model group: ".$self->id;

    return 1;
}

sub model_group_bridge_create {
    my $self = shift;

    my $model_group = $self->model_group;
    return 1 if not $model_group;
    return if not $model_group->uuid;

    my $project = Genome::Project->get($model_group->uuid);
    return 1 if not $project;

    my $model = $self->model;
    my $part = $project->parts(entity => $model);
    return 1 if $part;

    $part = $project->add_part(entity => $model);
    if ( not $part ) {
        die 'Failed to create project part for '.$project->id.' '.$model->id;
    }

    return 1;
}

sub model_group_bridge_delete {
    my $self = shift;
    $deleted_bridges{ $self->id }++;

    my $model_group = $self->model_group;
    return 1 if not $model_group;

    my $project = Genome::Project->get($model_group->uuid);
    return 1 if not $project;

    my $model = $self->model;
    my $part = $project->parts(entity => $model);
    return 1 if not $part;

    return 1 if $Genome::Site::TGI::Observers::Project::deleted_parts{ $part->id };
    if ( not $part->delete ) {
        die 'Failed to delete project part for '.$project->id.' '.$model->id;
    }

    return 1;
}

1;

