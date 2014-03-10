package Genome::Site::TGI::Synchronize::Classes::OrganismIndividual; 

use strict;
use warnings;

=pod
COMMON_NAME    VARCHAR2 (16)                    {null} {null}   ok
DESCRIPTION    VARCHAR2 (500)                   {null} {null}   ok
ETHNICITY      VARCHAR2 (64)                    {null} {null}   ok
FATHER_ID      NUMBER   (10)                    {null} {null}   ok
FULL_NAME      VARCHAR2 (64)                    {null} {null}   ok [name]
GENDER         VARCHAR2 (16)                    {null} {null}   ok
MOTHER_ID      NUMBER   (10)                    {null} {null}   ok
NAME           VARCHAR2 (64)                    {null} {null}   ok [upn]
NOMENCLATURE   VARCHAR2 (64)                    {null} {null}   ok
ORGANISM_ID    NUMBER   (10)                    {null} NOT NULL ok [id]
PARTICIPANT_ID VARCHAR2 (64)                    {null} {null}   NOT_SYNCED
RACE           VARCHAR2 (64)                    {null} {null}   ok
TAXON_ID       NUMBER   (10)                    {null} NOT NULL ok
13 properties, 12 copied, 11 updated
=cut

class Genome::Site::TGI::Synchronize::Classes::OrganismIndividual {
    is => 'Genome::Site::TGI::Synchronize::Classes::LimsBase',
    table_name => 'ORGANISM_INDIVIDUAL',
    id_by => [
        id => { is => 'Number', len => 10, column_name => 'ORGANISM_ID' },
    ],
    has => [
        name => { is => 'Text', column_name => 'FULL_NAME', }, # nullable in db
        taxon_id => { is => 'Number', },
    ],
    has_optional => [
        common_name => { is => 'Text', },
        description => { is => 'Text', },
        ethnicity => { is => 'Text', },
        father_id  => { is => 'Number', },
        gender => { is => 'Text', },
        mother_id  => { is => 'Number', },
        nomenclature => { is => 'Text', },
        race => { is => 'Text', },
        upn => { is => 'Text', column_name => 'NAME', },
    ],
    data_source => 'Genome::DataSource::Dwrac',
};

sub entity_name { return 'individual'; }

sub properties_to_copy {# 12
    return ( 'id', properties_to_keep_updated() );
}

sub properties_to_keep_updated {# 11
    return (qw/
        name
        common_name
        description
        ethnicity
        father_id 
        gender
        mother_id 
        nomenclature
        race
        taxon_id
        upn
    /);
}

sub lims_property_name_to_genome_property_name {
    my ($class, $name) = @_;
    my %lims_to_genome = (
        full_name => 'name',
        name => 'upn',
    );
    return $lims_to_genome{$name} if exists $lims_to_genome{$name};
    return $name;
}

1;

