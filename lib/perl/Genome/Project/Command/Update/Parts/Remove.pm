package Genome::Project::Command::Update::Parts::Remove;

use strict;
use warnings;

use Genome;

class Genome::Project::Command::Update::Parts::Remove { 
    is => 'Command::V2',
    has => [
        projects => {
            is => 'Genome::Project',
            is_many => 1,
            shell_args_position => 1,
            doc => 'Project(s), resolved via string.',
        },
        class_name => {
            is => 'Text',
            is_optional => 1,
            doc => 'The class name of the object to remove. To remove a value, leave undefined.',
        },
        value => {
            is => 'Text',
            doc => 'A string that can be any value or a value to get objects to remove from projects.',
        },
        role => {
            is => 'Text',
            is_optional => 1,
            doc => 'The role of the part',
        },
        label => {
            is => 'Text',
            is_optional => 1,
            doc => 'The label for the part',
        },
    ],
    doc => 'a part from projects',
};

sub help_synopsis {
    my $class = shift;
    my $command_name = $class->command_name;
    return <<HELP;
 Remove from a project (id 1) a project with id 2:
  $command_name id=1 --value id=2 --class-name Genome::Project
    
 Remove a priority value to project name 'High Priority':
  $command_name name='High Priority' --value 10 --label priority
HELP
}

sub help_detail {
    return 'Add values and objects to projects.'
}

sub execute {
    my $self = shift;

    $self->status_message('Remve from projects...');

    my @objects = ( not $self->class_name or $self->class_name eq 'UR::Value' )
    ? $self->_get_ur_value_for_value($self->value)
    : $self->_get_objects_for_class_and_value($self->class_name, $self->value);
    return if not @objects;

    my %params;
    for my $property (qw/ role label /) {
        my $value = $self->$property;
        $params{$property} = $value;
        $self->status_message(ucfirst($property).': '.( defined $value ? $value : 'NA' ));
    }

    for my $project ( $self->projects ) {
        $self->status_message('Project: '.$project->__display_name__);
        for my $object ( @objects ) {
            my $existing_part = $project->part(entity => $object, %params);
            if ( $existing_part ) {
                $existing_part->delete;
                $self->status_message('Removed: '.$object->__display_name__);
            }
            else {
                $self->status_message('Skipped, not assigned: '.$object->__display_name__);
            }
        }
    }

    $self->status_message('Done');

    return 1;
}

sub _get_objects_for_class_and_value {
    my ($self, $class, $value) = @_;

    $self->status_message('Class: '.$class);
    $self->status_message('Value: '.$class);

    my $bx = eval {
        UR::BoolExpr->resolve_for_string($class, $value);
    };
    if ( not $bx ) {
        $self->error_message($@);
        return;
    }

    my @objects = $class->get($bx);
    if ( not @objects ) {
        $self->status_message("No objects found for $class with value '$value'");
        return;
    }

    return @objects;
}

sub _get_ur_value_for_value {
    my ($self, $value) = @_;

    $self->status_message("Value for '$value'");

    if ( ref $value ) { # developer error
        $self->error_message('Cannot use a reference as value to get a UR::Value');
        return;
    }

    return UR::Value->get($value);
}

1;

