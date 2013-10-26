
package Genome::Site::TGI::Sample; 

# Adaptor for GSC Organism::Sample

# Do NOT use this module from anything in the GSC schema,
# though the converse will work just fine.

# This module should contain only UR class definitions,
# relationships, and support methods.

use strict;
use warnings;
use Genome;

=pod

    table_name => q|
        (
            select
                --fully precise and connected to LIMS
                s.organism_sample_id    id,
                s.full_name             name,
                s.common_name           common_name,

                -- collaborator's output
                s.sample_name           extraction_label,
                s.sample_type           extraction_type,
                s.description           extraction_desc,

                -- collaborator's source
                s.cell_type,
                s.tissue_label,
                s.tissue_name           tissue_desc,
                s.organ_name,

                -- patient, environment, or group for pools
                s.source_id,
                s.source_type,

		s.taxon_id

            from GSC.organism_sample s
        ) sample
    |,

=cut

class Genome::Site::TGI::Sample {
    is => 'Genome::Site::TGI::Measurable',
    table_name => 'ORGANISM_SAMPLE',
    id_by => [
        id                          => { is => 'Number',
                                        doc => 'the numeric ID for the specimen in both the LIMS and the analysis system', 
                                        column_name => 'ORGANISM_SAMPLE_ID',
                                    },
    ],
    has => [
        name                        => { is => 'Text',     len => 64, 
                                        doc => 'the fully qualified name for the sample (the "DNA NAME" in LIMS for both DNA and RNA)', 
                                        column_name => 'FULL_NAME',
                                    },
        subject_type => { is => 'Text', is_constant => 1, value => 'organism sample', column_name => '', },
        _nomenclature                => { column_name => 'NOMENCLATURE', default_value => "WUGC" }, 

    ],
    has_optional => [	
        taxon			    => { is => 'Genome::Site::TGI::Taxon', id_by => 'taxon_id' },
	default_genotype_seq_id     => {
            is => 'Text',
            doc => 'Seq ID of corresponding genotype data',
        },
        common_name                 => { is => 'Text', 
                                        doc => 'a name like "tumor1" for a given sample',                                        
                                    },

        extraction_label            => { is => 'Text', 
                                        doc => 'identifies the specimen sent from the laboratory which extracted DNA/RNA',
                                        column_name => 'SAMPLE_NAME',
                                    },
                
        extraction_type             => { is => 'Text', 
                                        doc => 'either "genomic dna" or "rna" in most cases', column_name => 'SAMPLE_TYPE' },
                
        extraction_desc             => { is => 'Text', 
                                        doc => 'notes specified when the specimen entered this site', 
                                        column_name => 'DESCRIPTION'
                                    },
                
        cell_type                   => { is => 'Text', len => 100,
                                        doc => 'typically "primary"' },

        tissue_label	            => { is => 'Text', 
                                        doc => 'identifies/labels the original tissue sample from which this extraction was made' },
        								
        tissue_desc                 => { is => 'Text', len => 64, 
                                        doc => 'describes the original tissue sample', column_name => 'TISSUE_NAME' },

        organ_name                  => { is => 'Text', len => 64, 
                                        doc => 'the name of the organ from which the sample was taken' }, 
        
        # these are optional only b/c our data is not fully back-filled
        source => { 
            is => 'Genome::Site::TGI::Measurable',
            id_by => 'source_id',
            where => [ 'subject_type in' => [qw/ organism_individual population_group /, 'organism individual', 'population group', ]],
            doc => 'The patient/individual organism from which the sample was taken, or the population for pooled samples.',
        },
        source_type                 => { is => 'Text',
                                        doc => 'either "organism individual" for individual patients, or "population group" for cross-individual samples' },
        
        source_name                 => { via => 'source', to => 'name' },
        
        source_common_name          => { via => 'source', to => 'common_name' },


        # the above are overly generic, since all of our sources are Genome::Individuals, and slow, so...
        patient                      => { is => 'Genome::Site::TGI::Individual', id_by => 'source_id',
                                           doc => 'The patient/individual organism from which the sample was taken.' },
        
        patient_name                 => { via => 'patient', to => 'name', doc => 'the system name for a patient (subset of the sample name)' },
        
        patient_common_name          => { via => 'patient', to => 'common_name', doc => 'names like AML1, BRC50, etc' },
        age => { 
            is => 'Number',
            via => 'attributes', 
            where => [ name => 'age', nomenclature => 'WUGC', ], 
            to => 'value',
            is_optional => 1,
            is_many => 0,
            is_mutable => 1,
            doc => 'Age of the patient at the time of sample taking.',
        },
        body_mass_index => {
            via => 'attributes',
            where => [ name => 'body_mass_index', nomenclature => 'WUGC', ] ,
            to => 'value',
            is_optional => 1,
            is_many => 0,
            is_mutable => 1,
            doc => 'BMI of the patient at the time of sample taking.',
        },
        tcga_name                   => { via => 'attributes', where => [ 'nomenclature like' => 'TCGA%', name => 'biospecimen_barcode_side'], to => 'value' },
        sub_type                    => { calculate_from => ['_sub_type1','_sub_type2'], calculate => q|$_sub_type1 or $_sub_type2| }, 
        _sub_type1                  => { via => 'attributes', where => [ name => 'sub-type' ], to => 'value' },
        _sub_type2                  => { via => 'attributes', where => [ name => 'subtype' ], to => 'value' },
        
        models                      => { is => 'Genome::Model', reverse_as => 'subject', is_many => 1 },
    ],
    has_many => [
        attributes                  => { is => 'Genome::Site::TGI::Sample::Attribute', reverse_as => 'sample', specify_by => 'name', is_optional => 1, is_many => 1, },
        libraries                   => { is => 'Genome::Library', reverse_id_by => 'sample' },
        solexa_lanes                => { is => 'Genome::InstrumentData::Solexa', reverse_id_by => 'sample' },
        solexa_lane_names           => { via => 'solexa_lanes', to => 'full_name' },
    ],
    doc         => 'a single specimen of DNA or RNA extracted from some tissue sample',
    data_source => 'Genome::DataSource::Dwrac',
};

sub __display_name__ {
    my $self = $_[0];
    return $self->name . ($self->patient_common_name ? ' (' . $self->patient_common_name . ' ' . $self->common_name . ')' : '');
}

sub create {
    my $class = shift;
    my $self = $class->SUPER::create(@_);
    return if not defined $self;

    # FIXME add tissue, nomenclature?

    return $self;
}

sub sample_type {
    shift->extraction_type(@_);
}

sub models {
    my $self = shift;
    my @m = Genome::Model->get(subject_id => $self->id, subject_class_name => $self->class);
    return @m;
}

1;

