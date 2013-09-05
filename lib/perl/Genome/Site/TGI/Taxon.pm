package Genome::Site::TGI::Taxon; 

use strict;
use warnings;

class Genome::Site::TGI::Taxon {
    is => 'Genome::Site::TGI::Measurable',
    table_name => "organism_taxon",
    id_by => [
        id         => { is => 'Number', len => 10, column_name => 'TAXON_ID' },
    ],
    has => [
        taxon_id                        => { calculate => q|$self->id| }, 
        name => { 
            is => "Text", 
            len => 99,
            column_name => 
            'SPECIES_NAME', 
            doc => 'Name (species name) of the taxon', 
        },
        domain => { is => "Text",   len => 9, valid_values => [qw/ Archea Bacteria Eukaryota Virus Unknown /], doc => 'Domain of the taxon.', },
        species_name                    => { is => "Text",   len => 64, 
                                                calculate => q|$name|, 
                                                calculate_from => ['name'],
                                                #TODO: this actually embeds the strain name, parse it away
                                            },
        subject_type => { is => 'Text', is_constant => 1, value => 'organism taxon', column_name => '', },
    ],
    has_optional => [
        strain_name                     => { is => "Text",   len => 32, },
        species_latin_name              => { is => "Text",   len => 64 },
        ncbi_taxon_id                   => { is => "Number", len => 10 },
        ncbi_taxon_species_name         => { is => "Text",   len => 128 },
        locus_tag                       => { is => "Text",   len => 200 },
        gram_stain_category             => { is => "Text",   len => 32, column_name => 'GRAM_STAIN_CATEGORY', valid_values => [qw/ positive negative indeterminate /], },
        estimated_genome_size           => { is => "Number", len => 12, column_name => 'ESTIMATED_ORGANISM_GENOME_SIZE' },
        current_default_org_prefix      => { is => "Text",   len => 2, is_mutable => 0 },
        current_genome_refseq_id        => { is => "Number", len => 15 },
        
        model_member                    => { is => 'Genome::Site::TGI::Individual', id_by => 'model_member_id',
                                            doc => 'the model individual or inbred group sequenced as a reference for this taxon' },
        
        model_member_id                 => { is => "Number", len => 10, column_name => 'MODEL_INDIVIDUAL_ORGANISM_ID' },
        
        _legacy_org_id                  => { is => "Number", len => 10, column_name => 'LEGACY_ORG_ID' },
        _next_amplicon_iteration         => { is => "Number", len => 8, column_name => 'NEXT_AMPLICON_ITERATION' },
    ],
    has_many_optional => [
        individuals                     => { is => 'Genome::Site::TGI::Individual', reverse_id_by => 'taxon',  
                                            doc => 'all tracked individual organisms (patients/research subjects) of this species/strain' },                         
 
        population_groups               => { is => 'Genome::Site::TGI::PopulationGroup', reverse_id_by => 'taxon',
                                            doc => 'all defined population groups for this species/strain' },

        members                         => {
                                            calculate => q|($self->individuals)|,
                                            doc => 'all individuals AND defined population groups' },


        samples                         => { is => 'Genome::Site::TGI::Sample', is_many => 1, reverse_id_by => 'taxon',
                                            # if we had complete tracking, and it were efficient, we'd get this via populations above
                                            doc =>  'all DNA/RNA extractions from associated individuals and population groups' },
    ],
    doc => 'a species, strain, or other taxonomic unit',
    data_source => 'Genome::DataSource::Dwrac',
};

sub create {
    Carp::confess('Can not create taxons in LIMS!');
}

1;

