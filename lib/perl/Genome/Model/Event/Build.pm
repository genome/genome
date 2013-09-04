package Genome::Model::Event::Build;

use strict;
use warnings;

use Genome;

# this stub exists only so command_name still returns the type_name expected by the database 

class Genome::Model::Event::Build {
    is => [ 'Genome::Model::Event' ],
    has => [
        data_directory => { via => 'build' },
    ],
};

sub command_name {
    my $class = ref($_[0]) || $_[0];
    return $class->SUPER::command_name unless $class eq __PACKAGE__;
    return 'genome model build';
}

sub command_name_brief {
    my $class = ref($_[0]) || $_[0];
    return $class->SUPER::command_name_brief unless $class eq __PACKAGE__;
    return 'build';
}

sub resolve_stage_name_for_class {
    my $self = shift;
    my $class = shift;
    my $pp = $self->model->processing_profile;
    for my $stage_name ($pp->stages) {
        my $found_class = grep { $class =~ /^$_/ } $pp->classes_for_stage($stage_name);
        if ($found_class) {
            return $stage_name;
        }
    }
    my $error_message = "No class found for '$class' in build ". $self->class ." stages:\n";
    for my $stage_name ($pp->stages) {
        $error_message .= $stage_name ."\n";
        for my $class ($pp->classes_for_stage($stage_name)) {
            $error_message .= "\t". $class ."\n";
        }
    }
    $self->error_message($error_message);
    return;
}

sub events_for_stage {
    my $self = shift;
    my $stage_name = shift;
    my $pp = $self->model->processing_profile;
    my @events;
    for my $class ($pp->classes_for_stage($stage_name)) {
        push @events, $self->events_for_class($class);
    }
    return @events;
}

sub events_for_class {
    my $self = shift;
    my $class = shift;

    my @class_events = $class->get(
                                   model_id => $self->model_id,
                                   build_id => $self->build_id,
                               );

    #Not sure if every class is supposed to have return of events
    #but adding the line below makes the tests pass for now
    return unless @class_events;

    return sort {$a->id cmp $b->id} @class_events;
}

1;
