package Genome::Observers::SampleLibrary;

use strict;
use warnings;

use Carp;

Genome::Sample->add_observer(
    aspect => 'delete',
    callback => \&delete_libraries
);

sub delete_libraries {
    my $self = shift;
    my @libraries = Genome::Library->get(sample_id => $self->id);
    if (@libraries) {
        Carp::confess "Cannot delete sample " . $self->__display_name__ . 
            ", there are " . scalar @libraries . " that are derived from it!";
    }
    return 1;
}
            
1;

