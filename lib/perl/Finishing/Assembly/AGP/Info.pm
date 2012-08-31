package Finishing::Assembly::AGP::Info;

use strict;
use warnings;

use base 'Finfo::Singleton';

use Data::Dumper;

sub attributes_for_component_type
{
    my ($self, $type) = @_;

    $self->_enforce_instance;
    
    return unless Finfo::Validate->validate
    (
        attr => 'component type',
        value => $type,
        isa => [ 'in_list', agp_component_types() ],
        msg => 'fatal',
    );

    return (qw/
        object object_begin object_end part_number component_type gap_length gap_type linkage
    /) if $type eq 'N' or $type eq 'U';

    return  (qw/ 
        object object_begin object_end part_number component_type
        component_id component_begin component_end orientation 
        /);
}

sub component_position
{
    return 4;
}

# agp components
sub agp_component_types
{
    return map { $_->[0] } agp_components_data();
}

sub agp_component_attributes
{
    return (qw/ component_type status agp_type /);
}

sub agp_components_data
{
    # component_type status agp_type
    return
    (
        [ 'A', 'Active Finishing', 'object', ],
        [ 'D', 'Draft HTG', 'object', ],
        [ 'F', 'Finished HTG', 'object', ],
        [ 'G', 'Whole Genome Finishing', 'object', ],
        [ 'N', 'Gap', 'gap', ],
        [ 'O', 'Other', 'object', ],
        [ 'P', 'Pre Draft', 'object', ],
        [ 'U', 'Gap with Unkown Size', 'gap', ],
        [ 'W', 'WGS Contig', 'object', ]
    );
}

# agp gap types
sub valid_agp_gap_types
{
    return map { $_->[0] } agp_gap_types_data();
}

sub agp_gap_type_attributes
{
    return 'type';
}

sub agp_gap_types_data
{
    return 
    ( 
        ['unsized' ],
        ['fragment' ],
        ['clone' ],
        ['contig' ],
        ['centromere' ],
        ['short_arm' ],
        ['heterochromatin' ],
        ['telomere' ],
        ['repeat' ],
    );
}

1;

#$HeadURL: svn+ssh://svn/srv/svn/gscpan/perl_modules/trunk/Finishing/Assembly/AGP/Info.pm $
#$Id: Info.pm 30518 2007-11-30 22:45:36Z ebelter $
