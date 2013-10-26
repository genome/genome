package Genome::Site::TGI::Synchronize::Classes::OrganismTaxon; 

use strict;
use warnings;

=pod
CURRENT_DEFAULT_ORG_PREFIX     VARCHAR2 (2)                     {null} {null}   ok
CURRENT_GENOME_REFSEQ_ID       NUMBER   (15)                    {null} {null}   ok
DOMAIN                         VARCHAR2 (9)                     {null} {null}   ok
ESTIMATED_ORGANISM_GENOME_SIZE NUMBER   (12)                    {null} {null}   ok [estimated_genome_size]
GRAM_STAIN_CATEGORY            VARCHAR2 (32)                    {null} {null}   ok
LEGACY_ORG_ID                  NUMBER   (10)                    {null} NOT NULL NOT_SYNCED
LOCUS_TAG                      VARCHAR2 (200)                   {null} {null}   ok
MODEL_INDIVIDUAL_ORGANISM_ID   NUMBER   (10)                    {null} {null}   ok [model_member_id]
NCBI_TAXON_ID                  NUMBER   (10)                    {null} {null}   ok
NCBI_TAXON_SPECIES_NAME        VARCHAR2 (128)                   {null} {null}   ok
NEXT_AMPLICON_ITERATION        NUMBER   (8)                     1      NOT NULL NOT_SYNCED
SPECIES_LATIN_NAME             VARCHAR2 (64)                    {null} {null}   ok
SPECIES_NAME                   VARCHAR2 (64)                    {null} NOT NULL ok [name]
STRAIN_NAME                    VARCHAR2 (32)                    {null} {null}   ok
TAXON_ID                       NUMBER   (10)                    {null} NOT NULL ok [id]
15 properties, 13 copied, 12 updated
=cut

class Genome::Site::TGI::Synchronize::Classes::OrganismTaxon {
    is => 'UR::Object',
    table_name => "organism_taxon",
    id_by => [
        id => { is => 'Number', column_name => 'TAXON_ID' },
    ],
    has => [
        name => { is => 'Text', column_name => 'SPECIES_NAME', },
    ],
    has_optional => [
        current_default_org_prefix => { is => "Text", },
        current_genome_refseq_id => { is => "Number", },
        domain => { is => 'Text', valid_values => [qw/ Archea Bacteria Eukaryota Virus Unknown /], },
        estimated_genome_size => { is => "Number", column_name => 'ESTIMATED_ORGANISM_GENOME_SIZE' },
        gram_stain_category => { is => "Text", valid_values => [qw/ positive negative indeterminate /], },
        locus_tag => { is => 'Text', },
        model_member_id => { is => "Number", column_name => 'MODEL_INDIVIDUAL_ORGANISM_ID' },
        ncbi_taxon_id => { is => "Number", },
        ncbi_taxon_species_name => { is => "Text", },
        species_latin_name => { is => "Text", },
        strain_name => { is => "Text", },
    ],
    data_source => 'Genome::DataSource::Dwrac',
};

sub properties_to_copy {# 13
    return ( 'id', properties_to_keep_updated() );
}

sub properties_to_keep_updated {# 12
    return (qw/ 
        name
        current_default_org_prefix
        current_genome_refseq_id
        domain
        estimated_genome_size
        gram_stain_category
        locus_tag
        model_member_id
        ncbi_taxon_id
        ncbi_taxon_species_name
        species_latin_name
        strain_name
     /);
}

sub lims_property_name_to_genome_property_name {
    my ($class, $name) = @_;
    my %lims_to_genome = (
        species_name => 'name',
        estimated_organism_genome_size => 'estimated_genome_size',
        model_individual_organism_id => 'model_member_id',
    );
    return $lims_to_genome{$name} if exists $lims_to_genome{$name};
    return $name;
}

1;

