package Genome::Project;

use strict;
use warnings;
use Genome;
use Class::ISA;

class Genome::Project {
    is => ['Genome::Notable','Genome::Searchable'],
    id_generator => '-uuid',
    id_by => [
        id => { is => 'Text', }
    ],
    has => [
        name => {
            is => 'Text',
            doc => 'Name of the project',
        },
        fixed_size_name => {
            is => 'Text',
            is_calculated => 1,
            calculate_from => ['name'],
            calculate => sub {
                my ($n) = @_;
                return $n if (length($n) <= 20);
                return substr($n, 0, 12) . '..' . substr($n, length($n) - 6);
            },
            doc => 'Name of the project, fixed in size for project box on webpages',
        },
        user_ids => {
            is => 'Genome::Sys::User',
            via => 'parts',
            to => 'entity_id',
            where => [ 'entity_class_name' => 'Genome::Sys::User' ],
            is_mutable => 0,
            is_many => 1,
        },
    ],
    has_many_optional => [
        parts => {
            is => 'Genome::ProjectPart',
            is_mutable => 1,
            reverse_as => 'project',
            doc => 'All the parts that compose this project',
        },
        creator => {
            is => 'Genome::Sys::User',
            is_calculated => 1,
            calculate_from => ['parts'],
            calculate => sub {
                my (@parts) = @_;
                for my $p (@parts) {
                    if ($p->role && $p->role eq 'creator') {
                        return $p->entity();
                    }
                }
                return undef;
            }
        },
        part_set => {
            is => 'Genome::ProjectPart::Set',
            is_calculated => 1,
        },
        parts_count => { 
            is => 'Number', 
            via => 'part_set', 
            to => 'count',
            doc => 'The number of parts associated with this project',
        },
        entities => {
            via => 'parts',
            to => 'entity',
            doc => 'All the objects to which the parts point',
        },
    ],
    table_name => 'GENOME_PROJECT',
    schema_name => 'GMSchema',
    data_source => 'Genome::DataSource::GMSchema',
    doc => 'A project, can contain any number of objects (of any type)!',
};

sub create {
    my $class = shift;
    
    my $self = eval { $class->SUPER::create(@_) };
    if ($@ or not $self) {
        my $msg = "Could not create new object of type $class!" .  ($@ ? " Reason: $@" : '');
        $class->error_message($msg);
        die $msg;
    }

    # Set creator
    my $user_name = Genome::Sys->username;
    my $creator = Genome::Sys::User->get(username => $user_name);
    if ( not $creator ) {
        $self->error_message("Found no user entry for $user_name, cannot create project. " .
            "Either create an entry for $user_name by running 'genome sys user create' or contact informatics.");
        $self->delete;
        return;
    }

    if ( not $self->add_part(entity => $creator, role => 'creator') ) {
        $self->error_message("Failed to add creater '$user_name' to project");
        $self->delete;
        return;
    }

    return $self;
}

sub get_part {
    my ($self, $obj) = @_;

    my @parts = Genome::ProjectPart->get(
        entity_class_name => $obj->class,
        entity_id => $obj->id
    );

    return @parts;
}

sub __display_name__ {
    my $self = shift;
    return $self->name."(".$self->id.")";
}

sub rename {
    my ($self, $new_name) = @_;

    unless ($new_name) {
        $self->error_message('No new name given to rename project');
        return;
    }

    my @projects = Genome::Project->get(name => $new_name);
    if (@projects) {
        $self->error_message("Failed to rename project (" . $self->id .
            ") from '" . $self->name . "' to '$new_name' because one already exists.");
        return;
    }

    my $old_name = $self->name;
    $self->name($new_name);
    my $rv = eval { $self->name($new_name) };
    if ($@ or not $rv) {
        $self->error_message("Could not rename project " . $self->__display_name__ .
            " from $old_name to $new_name!" .
            ($@ ? " Reason: $@" : ""));
        return;
    }

    $self->status_message("Renamed project from '$old_name' to '$new_name'");

    return 1;
}

sub get_parts_of_class {
    my $self = shift;
    my $desired_class = shift;
    croak $self->error_message('missing desired_class argument') unless $desired_class;

    my @parts = $self->parts;
    return unless @parts;

    my @desired_parts;
    for my $part (@parts) {
        my @classes = Class::ISA::self_and_super_path($part->entity->class);
        push @desired_parts, $part if grep { $_ eq $desired_class } @classes;
    }

    return @desired_parts;
}

sub delete {

    my ($self) = @_;

    for my $part ($self->parts) {
        $part->delete();
    }    

    return $self->SUPER::delete();
}



1;

