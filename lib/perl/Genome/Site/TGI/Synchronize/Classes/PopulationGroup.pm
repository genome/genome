package Genome::Site::TGI::Synchronize::Classes::PopulationGroup; 

use strict;
use warnings;

=pod
DESCRIPTION VARCHAR2 (500)                   {null} {null}   ok
NAME        VARCHAR2 (64)                    {null} {null}   ok
PG_ID       NUMBER   (10)                    {null} NOT NULL ok [id]
TAXON_ID    NUMBER   (10)                    {null} NOT NULL ok
4 properties, 4 copied, 3 updated
=cut

class Genome::Site::TGI::Synchronize::Classes::PopulationGroup { 
    is => 'UR::Object',
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
    data_source => 'Genome::DataSource::Dwrac',
};

sub properties_to_copy {# 4
    return ( 'id', properties_to_keep_updated() );
}

sub properties_to_keep_updated {# 3
    return (qw/
        name
        taxon_id
        description
    /);
}

sub lims_property_name_to_genome_property_name {
    return $_[1];
}

1;

