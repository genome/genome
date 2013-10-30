package Genome::Site::TGI::Synchronize::Classes::OrganismSample; 

use strict;
use warnings;

=pod
ARRIVAL_DATE                DATE                             {null} {null}   NOT_SYNCED
CELL_TYPE                   VARCHAR2 (100)                   {null} NOT NULL ok
COMMON_NAME                 VARCHAR2 (32)                    {null} {null}   ok
CONFIRM_EXT_GENOTYPE_SEQ_ID NUMBER   (15)                    {null} {null}   NOT_SYNCED
DEFAULT_GENOTYPE_SEQ_ID     NUMBER   (15)                    {null} {null}   ok [updated by reconcile genotype]
DESCRIPTION                 VARCHAR2 (500)                   {null} {null}   ok [extraction_desc]
DNA_PROVIDER_ID             NUMBER   (10)                    {null} {null}   NOT_SYNCED
FULL_NAME                   VARCHAR2 (64)                    {null} {null}   ok [name]
GENDER                      VARCHAR2 (16)                    {null} {null}   NOT_SYNCED
GENERAL_RESEARCH_CONSENT    NUMBER   (1)                     {null} {null}   NOT_SYNCED
IS_PROTECTED_ACCESS         NUMBER   (1)                     {null} {null}   NOT_SYNCED
IS_READY_FOR_ANALYSIS       NUMBER   (1)                     {null} {null}   NOT_SYNCED
NOMENCLATURE                VARCHAR2 (64)                    {null} NOT NULL NOT_SYNCED
ORGANISM_SAMPLE_ID          NUMBER   (20)                    {null} NOT NULL ok [id]
ORGAN_NAME                  VARCHAR2 (64)                    {null} {null}   ok
REFERENCE_SEQUENCE_SET_ID   NUMBER   (10)                    {null} {null}   NOT_SYNCED
SAMPLE_NAME                 VARCHAR2 (64)                    {null} {null}   ok [extraction_label]
SAMPLE_TYPE                 VARCHAR2 (32)                    {null} {null}   ok [extraction_type]
SOURCE_ID                   NUMBER   (10)                    {null} {null}   ok
SOURCE_TYPE                 VARCHAR2 (64)                    {null} {null}   NOT_SYNCED
TAXON_ID                    NUMBER   (10)                    {null} {null}   NOT_SYNCED
TISSUE_LABEL                VARCHAR2 (64)                    {null} {null}   ok
TISSUE_NAME                 VARCHAR2 (64)                    {null} {null}   ok [tissue_desc]
24 properties, 15 to copy, 12 to update (exclude id, name, default_genotype_seq_id)
=cut
#TODO SAMPLE_ATTRIBUTES

class Genome::Site::TGI::Synchronize::Classes::OrganismSample {
    is => 'Genome::Site::TGI::Synchronize::Classes::LimsBase',
    table_name => 'ORGANISM_SAMPLE',
    id_by => [
        id => { is => 'Number', column_name => 'ORGANISM_SAMPLE_ID', },
    ],
    has => [
        cell_type => { is => 'Text', },
        name => { is => 'Text', column_name => 'FULL_NAME', },
        nomenclature => { is => 'Text', }, 
    ],
    has_optional => [	
        common_name             => { is => 'Text', },
        default_genotype_seq_id => { is => 'Text', },
        extraction_desc         => { is => 'Text', column_name => 'DESCRIPTION', },
        extraction_label        => { is => 'Text', column_name => 'SAMPLE_NAME', },
        extraction_type         => { is => 'Text', column_name => 'SAMPLE_TYPE', },
        organ_name              => { is => 'Text', },
        source_id               => { is => 'Number', },
        tissue_desc             => { is => 'Text', column_name => 'TISSUE_NAME', },
        tissue_label	        => { is => 'Text', },
    ],
    data_source => 'Genome::DataSource::Dwrac',
};

sub entity_name { return 'sample'; }

sub properties_to_copy {# 15
    return ( 'id', 'name', 'default_genotype_seq_id', properties_to_keep_updated() );
}

sub properties_to_keep_updated {# 12
    return (qw/ 
        cell_type
        common_name
        extraction_desc
        extraction_label
        extraction_type
        organ_name
        source_id
        tissue_desc
        tissue_label
        /);
}

sub lims_property_name_to_genome_property_name {
    my ($class, $name) = @_;
    my %lims_to_genome = (
        full_name => 'name',
        description => 'extraction_desc',
        sample_name => 'extraction_label',
        sample_type => 'extraction_type',
        tissue_name => 'tissue_desc',
    );
    return $lims_to_genome{$name} if exists $lims_to_genome{$name};
    return $name;
}

1;

