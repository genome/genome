package Genome::Site::TGI::Synchronize::Classes::PopulationGroup; 

use strict;
use warnings;

=pod
DESCRIPTION VARCHAR2 (500)                   {null} {null}   ok
NAME        VARCHAR2 (64)                    {null} {null}   ok
PG_ID       NUMBER   (10)                    {null} NOT NULL ok [id]
TAXON_ID    NUMBER   (10)                    {null} NOT NULL ok
4 properties, 4 copied, 3 updated
Example: 2774516127
=cut

class Genome::Site::TGI::Synchronize::Classes::PopulationGroup { 
    is => 'Genome::Site::TGI::Synchronize::Classes::LimsBase',
    table_name => 'POPULATION_GROUP',
    id_by => [
        id => { is => 'Number', column_name => 'PG_ID' },
    ],
    has => [
        name => { is => 'Text', }, # nullable in db
        taxon_id => { is => 'Number', },
    ],
    has_optional => [
        description => { is => 'Text', },
    ],
    has_many_optional => [
        member_ids => {
            is => 'Number',
            to => 'member_id',
            via => 'member_associations',
        },
        member_associations => { 
            is => 'Genome::Site::TGI::Synchronize::Classes::PopulationGroupMember', 
            reverse_id_by => 'population_group',
        },
    ],
    data_source => 'Genome::DataSource::Dwrac',
};

sub entity_name { return 'population group'; }

sub properties_to_copy {# 4
    return ( 'id', 'member_ids', properties_to_keep_updated() );
}

sub properties_to_keep_updated {# 3
    return (qw/
        name
        taxon_id
        description
    /);
}

1;

