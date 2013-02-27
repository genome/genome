package Genome::Site::TGI::Synchronize::Classes::MiscUpdate::RegionIndex454; 

use strict;
use warnings;

class Genome::Site::TGI::Synchronize::Classes::MiscUpdate::RegionIndex454 {
    is => 'Genome::Site::TGI::Synchronize::Classes::MiscUpdate::BaseForAttributes',
};

sub _after_update_value {
    my $self = shift;

    my $genome_entity = $self->genome_entity;
    return if not $genome_entity;

    my $genome_property_name = $self->genome_property_name;
    return if not $genome_property_name;

    return 1 if not grep { $genome_property_name eq $_ } (qw/ region_number index_sequence /);

    my $region_number = $genome_entity->region_number;
    return 1 if not defined $region_number; # should not happen

    my $index_sequence = $genome_entity->index_sequence;
    my $new_subset_name = ( defined $index_sequence )
    ? $region_number.'-'.$index_sequence
    : $region_number;
    $genome_entity->subset_name($new_subset_name);

    return 1;
}

1;

