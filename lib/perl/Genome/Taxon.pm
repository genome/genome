package Genome::Taxon;

use strict;
use warnings;
use Genome;

class Genome::Taxon {
    is => ['Genome::Subject','Genome::Searchable'],
    has => [
        taxon_id => {
            calculate => q|$self->id|
        },
        domain => {
            is => 'Text',
            via => 'attributes',
            to => 'attribute_value',
            where => [ attribute_label => 'domain' ],
            is_mutable => 1,
            valid_values => [qw/ Archaea Bacteria Eukaryota Unknown Virus /],
            doc => 'Domain of this taxon, eg eukaryota',
        },
        species_name => {
            is => 'Text',
            calculate_from => 'name',
            calculate => q{
                $self = shift;
                my $new_name = shift;
                if ($new_name) { $self->name($new_name); }
                return $name;
            },
            doc => 'Plain text species name',
        },
        subject_type => {
            is_constant => 1,
            is_classwide => 1,
            value => 'species_name',
        },
    ],
    has_optional => [
#       This was once on the Subject base class but this caused UR to not know how to join over to
#       SubjectAttributes.
        common_name => {
            calculate_from => 'name',
            calculate => q{
                $self = shift;
                my $new_name = shift;
                if ($new_name) { $self->name($new_name); }
                return $name
            },
        },
        strain_name=> {
            is => 'Text',
            via => 'attributes',
            to => 'attribute_value',
            where => [ attribute_label => 'strain_name' ],
            is_mutable => 1,
            doc => 'Name of the strain, if applicable',
        },
        species_latin_name => {
            is => 'Text',
            via => 'attributes',
            to => 'attribute_value',
            where => [ attribute_label => 'species_latin_name' ],
            is_mutable => 1,
            doc => 'Latin species name',
        },
        ncbi_taxon_id  => {
            is => 'Number',
            via => 'attributes',
            to => 'attribute_value',
            where => [ attribute_label => 'ncbi_taxon_id' ],
            is_mutable => 1,
            doc => 'NCBI taxon ID, if available',
        },
        ncbi_taxon_species_name => {
            is => 'Text',
            via => 'attributes',
            to => 'attribute_value',
            where => [ attribute_label => 'ncbi_taxon_species_name' ],
            is_mutable => 1,
            doc => 'NCBI taxon species name, if available',
        },
        locus_tag => {
            is => 'Text',
            via => 'attributes',
            to => 'attribute_value',
            where => [ attribute_label => 'locus_tag' ],
            is_mutable => 1,
            doc => 'Locus tag of the taxon',
        },
        gram_stain_category => {
            is => 'Text',
            via => 'attributes',
            to => 'attribute_value',
            where => [ attribute_label => 'gram_stain_category' ],
            is_mutable => 1,
            valid_values => [qw/ indeterminate variable positive negative /],
            doc => 'Gram stain of the taxon, if applicable',
        },
        estimated_genome_size => {
            is => 'Number',
            via => 'attributes',
            to => 'attribute_value',
            where => [ attribute_label => 'estimated_genome_size' ],
            is_mutable => 1,
            doc => 'Estimated size of this taxon\'s genome',
        },
        current_default_org_prefix => {
            is => 'Text',
            via => 'attributes',
            to => 'attribute_value',
            where => [ attribute_label => 'current_default_org_prefix' ],
            is_mutable => 1,
            doc => 'Current organism prefix',
        },
        current_genome_refseq_id => {
            is => 'Number',
            via => 'attributes',
            to => 'attribute_value',
            where => [ attribute_label => 'current_genome_refseq_id' ],
            is_mutable => 1,
        },
        # TODO Should remove these, used by XML view though
        model_member => {
            is => 'Genome::Individual',
            id_by => 'model_member_id',
            doc => 'the model individual or inbred group sequenced as a reference for this taxon'
        },
        model_member_id => {
            is => 'Number',
            via => 'attributes',
            to => 'attribute_value',
            where => [ attribute_label => 'model_member_id' ],
            is_mutable => 1,
        },
        _legacy_org_id => {
            is => 'Number',
            via => 'attributes',
            to => 'attribute_value',
            where => [ attribute_label => '_legacy_org_id' ],
            is_mutable => 1,
        },
        _next_amplicon_iteration => {
            is => 'Number',
            via => 'attributes',
            to => 'attribute_value',
            where => [ attribute_label => '_next_amplicon_iteration' ],
            is_mutable => 1,
        },
    ],
    has_many_optional => [
        individuals => {
            is => 'Genome::Individual',
            reverse_id_by => 'taxon',
            doc => 'All tracked individual organisms (patients/research subjects) of this species/strain'
        },
        population_groups => {
            is => 'Genome::PopulationGroup',
            reverse_id_by => 'taxon',
            doc => 'All defined population groups for this species/strain'
        },
        members => {
            calculate => q|($self->individuals)|,
            doc => 'All individuals AND defined population groups'
        },
        samples => {
            is => 'Genome::Sample',
            reverse_id_by => 'taxon',
            doc => 'All DNA/RNA extractions from associated individuals and population groups' },
    ],
    doc => 'A species, strain, or other taxonomic unit',
};

sub __display_name__ {
    my $self = $_[0];
    return $self->name . ' (' . $self->id . ')';
}

sub get_source {
    return;
}

1;

