package Genome::Sample;

use strict;
use warnings;

use Genome;

my $default_nomenclature = 'GMS'; 

class Genome::Sample {
    is => ['Genome::Subject','Genome::Searchable'],
    has => [
        sample_id => {
            is => 'Text',
            calculate_from => 'subject_id',
            calculate => q{ return $subject_id },
        },
        subject_type => {
            is_constant => 1,
            is_classwide => 1,
            value => 'sample_name'
        },
    ],
    has_optional => [
        common_name => {
            is => 'Text',
            via => 'attributes',
            to => 'attribute_value',
            where => [ attribute_label => 'common_name', nomenclature => $default_nomenclature ],
            is_mutable => 1,
            doc => 'Typically tumor, normal, etc. A very brief description of the sample',
        },
        individual_common_name => {
            is => 'Text',
            via => 'source',
            to => 'common_name',
            doc => 'AML45, BRC1, etc'
        },
        extraction_label => {
            is => 'Text',
            via => 'attributes',
            to => 'attribute_value',
            where => [ attribute_label => 'extraction_label', nomenclature => $default_nomenclature ],
            is_mutable => 1,
            doc => 'Identifies the specimen sent from the laboratory which extracted DNA/RNA',
        },
        extraction_type => {
            is => 'Text',
            via => 'attributes',
            to => 'attribute_value',
            where => [ attribute_label => 'extraction_type', nomenclature => $default_nomenclature ],
            is_mutable => 1,
            doc => 'Either "genomic dna" or "rna" in most cases',
        },
        sample_type => {
            calculate_from => 'extraction_type',
            calculate => q{
                $self = shift;
                my $new_type = shift;
                if ($new_type) { $self->extraction_type($new_type); }
                return $extraction_type;
            },
        },
        is_rna => {
            calculate_from => [qw/ extraction_type /],
            calculate => q|
                return if not defined $extraction_type;
                for my $rna_sample_type ( 'rna', 'cdna', 'total rna', 'cdna library', 'mrna', 'pooled rna' ) {
                    return 1 if $extraction_type eq $rna_sample_type;
                }
                return;
            |,
        },
        extraction_desc => {
            is => 'Text',
            via => 'attributes',
            to => 'attribute_value',
            where => [ attribute_label => 'extraction_desc', nomenclature => $default_nomenclature ],
            is_mutable => 1,
            doc => 'Notes specified when the specimen entered this site',
        },
        cell_type => {
            is => 'Text',
            via => 'attributes',
            to => 'attribute_value',
            where => [ attribute_label => 'cell_type', nomenclature => $default_nomenclature ],
            is_mutable => 1,
            doc => 'Typically "primary"'
        },
        tissue_label => {
            is => 'Text',
            via => 'attributes',
            to => 'attribute_value',
            where => [ attribute_label => 'tissue_label', nomenclature => $default_nomenclature ],
            is_mutable => 1,
            doc => 'Identifies/labels the original tissue sample from which this extraction was made'
        },
        tissue_desc => {
            is => 'Text',
            via => 'attributes',
            to => 'attribute_value',
            where => [ attribute_label => 'tissue_desc', nomenclature => $default_nomenclature ],
            is_mutable => 1,
            doc => 'Describes the original tissue sample',
        },
        organ_name => {
            is => 'Text',
            via => 'attributes',
            to => 'attribute_value',
            where => [ attribute_label => 'organ_name', nomenclature => $default_nomenclature ],
            is_mutable => 1,
            doc => 'The name of the organ from which the sample was taken'
        },
        disease => {
            is => 'Text',
            via => 'attributes',
            to => 'attribute_value',
            where => [ attribute_label => 'disease', nomenclature => $default_nomenclature ],
            is_mutable => 1,
            doc => 'The name of the disease if present in the sample.',
        },
        default_genotype_data_id => {
            is => 'Number',
            via => 'attributes',
            to => 'attribute_value',
            where => [ attribute_label => 'default_genotype_data', nomenclature => $default_nomenclature ],
            is_mutable => 1,
            doc => 'ID of genotype microarray data associated with this sample',
        },
        default_genotype_data => {
            is => 'Genome::InstrumentData::Imported',
            id_by => 'default_genotype_data_id',
            doc => 'Genotype microarray instrument data object',
        },
        source_id => {
            is => 'Number',
            via => 'attributes',
            to => 'attribute_value',
            where => [ attribute_label => 'source_id', nomenclature => $default_nomenclature ],
            is_mutable => 1,
            doc => 'ID of the source of this sample, either a Genome::Individual or Genome::PopulationGroup',
        },
        source => {
            is => 'Genome::SampleSource',
            id_by => 'source_id',
            doc => 'The patient/individual organism or group from which the sample was taken, or the population for pooled samples.',
        },
        source_type => {
            is => 'Text',
            calculate_from => 'source',
            calculate => q{
                return unless $source;
                return $source->subject_type;
            },
            doc => 'Plain text type of the sample source',
        },
        source_name => {
            via => 'source',
            to => 'name',
            doc => 'Name of the sample source',
        },
        source_common_name => {
            via => 'source',
            to => 'common_name',
            doc => 'Common name of the sample source',
        },
        # These patient properties are for convenience, since the vast majority of sample
        # sources are of type Genome::Individual
        patient => {
            is => 'Genome::Individual',
            id_by => 'source_id',
            doc => 'The patient/individual organism from which the sample was taken.'
        },
        patient_name => {
            via => 'patient',
            to => 'name',
            doc => 'The system name for a patient (subset of the sample name)'
        },
        patient_common_name => {
            via => 'patient',
            to => 'common_name',
            doc => 'Common name of the patient, eg AML1',
        },
        age => {
            is => 'Number',
            via => 'attributes',
            to => 'attribute_value',
            where => [ attribute_label => 'age' ],
            is_mutable => 1,
            doc => 'Age of the patient at the time of sample taking.',
        },
        body_mass_index => {
            is => 'Text',
            via => 'attributes',
            to => 'attribute_value',
            where => [ attribute_label => 'body_mass_index' ],
            is_mutable => 1,
            doc => 'BMI of the patient at the time of sample taking.',
        },
        tcga_name => {
            via => 'attributes',
            to => 'attribute_value',
            where => [ 'nomenclature like' => 'TCGA%', attribute_label => 'biospecimen_barcode_side'],
            is_mutable => 1,
            doc => 'TCGA name of the sample, if available',
        },
        taxon => {
            is => 'Genome::Taxon',
            via => 'source',
            to => 'taxon',
            doc => 'Taxon for this sample via the source.',
        },
        taxon_id => {
            is => 'Number',
            via => 'source',
            to => 'taxon_id',
            doc => 'Taxon id for this sample via the source.',
        },
        species_name => {
            via => 'taxon',
            to =>  'name',
            doc => 'Name of the species of the sample source\'s taxonomic category'
        },
        sub_type => {
            calculate_from => ['_sub_type1','_sub_type2'],
            calculate => q|$_sub_type1 or $_sub_type2|
        },
        _sub_type1 => {
            via => 'attributes',
            to => 'attribute_value',
            where => [ attribute_label => 'sub-type' ],
            is_mutable => 1,
        },
        _sub_type2 => {
            via => 'attributes',
            to => 'attribute_value',
            where => [ attribute_label => 'subtype' ],
            is_mutable => 1,
        },
    ],
    has_many_optional => [
        models => {
            is => 'Genome::Model',
            reverse_as => 'subject',
            doc => 'Models that use this sample',
        },
        libraries => {
            is => 'Genome::Library',
            reverse_as => 'sample',
            is_many => 1,
            doc => 'Libraries that were created from the sample',
        },
        instrument_data => {
            is => 'Genome::InstrumentData',
            via => 'libraries',
            is_many => 1,
            doc => 'Instrument data from all DNA libraries from this sample',
        },
    ],
    doc => 'A single specimen of DNA or RNA extracted from some tissue sample',
};

sub get_source {
    my $self = shift;
    return $self->source;
}

sub __display_name__ {
    my $self = $_[0];
    return $self->name . ' (' . (($self->source && $self->source->common_name) ? $self->source->common_name . ($self->common_name ? ' ' . $self->common_name  : '') . ' ' : '') . $self->id . ')';
}

sub check_genotype_data {
    my $self = shift;
    my $genotype_instrument_data = shift;

    Carp::confess $self->error_message("No genotype instrument data provided.")
        unless $genotype_instrument_data;

    Carp::confess $self->error_message("Genotype instrument data is not a Genome::InstrumentData::Imported object.")
        unless $genotype_instrument_data->isa('Genome::InstrumentData::Imported');

    Carp::confess $self->error_message("Instrument data is not a 'genotype file' format.")
       unless $genotype_instrument_data->import_format && $genotype_instrument_data->import_format eq 'genotype file';

    my $genotype_sample = $genotype_instrument_data->sample;
    my $genotype_source = $genotype_sample->source;
    unless ($self->source->class eq $genotype_source->class and $self->source->id eq $genotype_source->id) {
        Carp::confess $self->error_message("Genotype instrument data has source " . $genotype_source->__display_name__ .
            " but sample has source " . $self->source->__display_name__);
    }

    return 1;
}

sub set_default_genotype_data {
    my ($self, $genotype_data_id, $allow_overwrite) = @_;
    $allow_overwrite ||= 0;
    Carp::confess 'Not given genotype instrument data to assign to sample ' . $self->id unless $genotype_data_id;

    unless ($genotype_data_id eq 'none') {
        my $genotype_instrument_data = Genome::InstrumentData::Imported->get($genotype_data_id);
        Carp::confess "Could not find any instrument data with id $genotype_data_id!" unless $genotype_instrument_data;
        Carp::confess "Genotype instrument data $genotype_data_id is not valid!"
            unless $self->check_genotype_data($genotype_instrument_data);
    }

    if (defined $self->default_genotype_data_id) {
        return 1 if $self->default_genotype_data_id eq $genotype_data_id;
        unless ($allow_overwrite) {
            Carp::confess "Attempted to overwrite current genotype instrument data id " . $self->default_genotype_data_id .
                " for sample " . $self->id . " with genotype data id $genotype_data_id " .
                " without setting the overwrite flag!";
        }
    }

    $self->default_genotype_data_id($genotype_data_id);

    my @genotype_models = $self->default_genotype_models;
    unless (@genotype_models) {
        $self->warning_message("Found no default genotype models using sample " . $self->__display_name__);
    }
    else {
        for my $genotype_model ($self->default_genotype_models) {
            $genotype_model->request_builds_for_dependent_cron_ref_align;
            $self->status_message("Requested builds for reference alignment models dependent on genotype model " . $genotype_model->id);
        }
    }

    return 1;
}

sub default_genotype_models {
    my $self = shift;

    my $genotype_data_id = $self->default_genotype_data_id;
    return unless defined $genotype_data_id;
    return if $genotype_data_id eq 'none';

    my $genotype_data = $self->default_genotype_data;
    return unless $genotype_data;

    my @inputs = Genome::Model::Input->get(
        value_class_name => $genotype_data->class,
        value_id => $genotype_data->id,
        name => 'instrument_data',
    );
    my @models = map { $_->model } @inputs;
    @models = grep { $_->subclass_name eq 'Genome::Model::GenotypeMicroarray' } @models;

    return @models;
}

sub delete {
    my $self = shift;

    for my $library ( $self->libraries ) {
        $library->delete;
    }

    return $self->SUPER::delete;
}

1;

