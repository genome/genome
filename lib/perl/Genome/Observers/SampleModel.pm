package Genome::Observers::SampleModel;

use strict;
use warnings;
use Carp;

UR::Observer->register_callback(
    subject_class_name => 'Genome::Sample',
    aspect => 'delete',
    callback => \&delete_models,
);

sub delete_models {
    my $self = shift;
    my @models = Genome::Model->get(subject_id => $self->id);
    if (@models) {
        Carp::confess "Cannot delete sample " . $self->__display_name__ . 
            ", there are " . scalar @models . " that have it as a subject!";
    }
    return 1;
}

1;

