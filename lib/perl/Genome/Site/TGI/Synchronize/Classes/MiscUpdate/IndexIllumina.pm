package Genome::Site::TGI::Synchronize::Classes::MiscUpdate::IndexIllumina; 

use strict;
use warnings;

class Genome::Site::TGI::Synchronize::Classes::MiscUpdate::IndexIllumina {
    is => 'Genome::Site::TGI::Synchronize::Classes::MiscUpdate::BaseForAttributes',
};

sub _after_update_value {
    my $self = shift;

    my $genome_property_name = $self->genome_property_name;
    return 1 if not grep { $genome_property_name eq $_ } (qw/  lane index_sequence /);

    my $genome_entity = $self->genome_entity;
    my $lane = $genome_entity->lane;
    return 1 if not defined $lane;

    my $index_sequence = $genome_entity->index_sequence;
    my $new_subset_name = ( defined $index_sequence )
    ? $lane.'-'.$index_sequence
    : $lane;
    $genome_entity->subset_name($new_subset_name);

    return 1;
}

1;

