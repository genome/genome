package Genome::Site::TGI::Observers::Project;

use strict;
use warnings;

our %groups_with_deleted_projects;
our %deleted_parts;

Genome::Project->add_observer(
    aspect => 'create',
    callback => \&create_callback,
);

Genome::Project->add_observer(
    aspect => 'delete',
    callback => \&delete_callback,
);

Genome::Project->add_observer(
    aspect => 'name',
    callback => \&project_rename,
);

Genome::ProjectPart->add_observer(
    aspect => 'create',
    callback => \&project_part_create,
);

Genome::ProjectPart->add_observer(
    aspect => 'delete',
    callback => \&project_part_delete,
);

sub create_callback {
    my ($self, $construction_method) = @_;
    my $class = $self->class;

    # Get project creator
    my $user_name = Genome::Sys->username;
    my $creator = Genome::Sys::User->get(username => $user_name);
    unless ($creator) {
        $self->delete;
        die "Failed to create project, could not find user $user_name";
    } # DO NOT SET HERE...IT IS DONE IN PROJECT CREATE

    # Make sure name is unique. Fail if name is not unique, but rename any conflicting
    # projects if this project is being created by apipe-builder.
    my ($existing_project) = grep { $self->id ne $_->id } $class->get(name => $self->name);
    if ($existing_project) {
        my $existing_project_creator_part = $existing_project->parts(role => 'creator');
        my $existing_project_creator_username = $existing_project_creator_part->entity->username;
        if ($user_name ne 'apipe-builder' or $existing_project_creator_username eq 'apipe-builder' ) {
            my $name = $self->name;
            $self->delete;
            die "There is already a project name '".$existing_project->name."', created by " . 
                $existing_project_creator_username . ". Select another name";
        }

        my $i = 0;
        my $old_name = $existing_project->name;
        my $new_name;
        do { 
            $new_name = $existing_project_creator_username. ' ' . $old_name .
            ($i ? '-'.$i : '');
            $i++;
        } while $class->get(name => $new_name);

        $self->status_message("There is another project with name " . $self->name . 
            " created by " . $existing_project_creator_username . 
            ", it will be renamed to $new_name");
        $existing_project->rename($new_name);
        # TODO email user about name change?
    }

    # Make sure corresponding model group exists
    my ($model_group) = Genome::ModelGroup->get(uuid => $self->id);
    if (not $model_group) {
        $self->status_message('Creating associated model group');
        $model_group = Genome::ModelGroup->create(
            name => $self->name,
            user_name => $creator->email,
            uuid => $self->id,
        );
        unless ($model_group) {
            $self->delete;
            die "Failed to create corresponding model group for project!";
        }
        $self->status_message('Create corresponding model group: '.$model_group->id);
    }

    if ( $self->name ne $model_group->name ) {
        die 'Project and model group names are not the same! '.join(' ', 'Project: ', map({ $self->$_ } (qw/ id name /)), 'Model Group: ', map({ $model_group->$_ } (qw/ id name /)));
    }

    return 1;
}

sub project_rename {
   my ($self, $property_name, $old_name, $new_name) = @_;

    my $model_group = Genome::ModelGroup->get(uuid => $self->id);
    return 1 if not $model_group;

    return 1 if $model_group->name eq $new_name;

    $self->status_message('Rename associated model group');
    $model_group->rename($new_name);

    return 1;
}

sub delete_callback {
    my $self = shift;

    $groups_with_deleted_projects{$self->id}++;
    return 1 if $Genome::Site::TGI::Observers::ModelGroup::projects_with_deleted_model_groups{$self->id};

    my ($model_group) = Genome::ModelGroup->get(uuid => $self->id);
    return 1 if not $model_group;

    $self->status_message('Deleting associated model group: '.$model_group->id);
    $model_group->delete;

    return 1;
}

sub project_part_create {
    my $self = shift;

    return 1 unless $self->entity_class_name->isa('Genome::Model');

    my ($model_group) = Genome::ModelGroup->get(uuid => $self->project_id);
    return 1 if not $model_group;

    my $model = $self->entity;
    my $bridge = $model_group->model_bridges(model => $model);
    return 1 if $bridge;

    $bridge = $model_group->add_model_bridge(model => $model);
    if ( not $bridge ) {
        die 'Failed to create model group bridge for '.$model_group->id.' '.$model->id;
    }

    return 1;
}

sub project_part_delete {
    my $self = shift;

    return 1 unless $self->entity_class_name->isa('Genome::Model');

    $deleted_parts{ $self->id }++;

    my $project = $self->project;
    my ($model_group) = Genome::ModelGroup->get(uuid => $self->project_id);
    return 1 if not $model_group;

    my $model = $self->entity;
    my $bridge = $model_group->model_bridges(model => $model);
    return 1 if not $bridge;
    return 1 if $Genome::Site::TGI::Observers::ModelGroup::deleted_bridges{ $bridge->id };

    if ( not $bridge->delete ) {
        die 'Failed to delete model group bridge for '.$model_group->id.' '.$model->id;
    }

    return 1;
}

1;

