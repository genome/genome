package Genome::InstrumentData;

use strict;
use warnings;

use Genome;

class Genome::InstrumentData {
    is => 'Genome::Notable',
    table_name => 'instrument.data',
    is_abstract => 1,
    subclass_description_preprocessor => __PACKAGE__ . '::_preprocess_subclass_description',
    attributes_have => [
        is_attribute => {
            is => 'Boolean',
            is_optional => 1,
        },
    ],
    subclassify_by => 'subclass_name',
    id_generator => '-uuid',
    id_by => [
        id => { is => 'Text', len => 64, },
    ],
    has => [
        seq_id => {
            calculate_from => 'id',
            calculate => q( return $id ),
            is_deprecated => 1,
        },
        subclass_name => {
            is => 'Text',
            len => 64,
        },
        sequencing_platform => {
            is => 'Text',
            len => 64,
        },
        import_format => {
            is => 'Text',
            via => 'attributes',
            to => 'attribute_value',
            is_mutable => 1,
            where => [ attribute_label => 'import_format' ],
        },
        import_source_name => {
            is => 'Text',
            via => 'attributes',
            to => 'attribute_value',
            is_mutable => 1,
            where => [ attribute_label => 'import_source_name' ],
        },
        description => {
            is => 'Text',
            via => 'attributes',
            to => 'attribute_value',
            is_mutable => 1,
            where => [ attribute_label => 'description' ],
        },
        library => {
            is => 'Genome::Library',
            id_by => 'library_id',
        },
        library_name => {
            via => 'library',
            to => 'name',
            is_deprecated => 1,
        },
        sample_id => {
            is => 'Text',
            via => 'library',
        },
        sample => {
            is => 'Genome::Sample',
            id_by => 'sample_id',
        },
        individual => {
            is => 'Genome::Individual',
            via => 'sample',
            to => 'patient',
        },
        taxon => {
            is => 'Genome::Taxon',
            via => 'sample'
        },
        sample_name => {
            via => 'sample',
            to => 'name',
            is_deprecated => 1,
        },
        models => {
            is => 'Genome::Model',
            reverse_as => 'instrument_data',
            is_many => 1,
        },
        deprecated_source_name => {
            is => 'Text',
            len => 64,
            is_optional => 1,
            column_name => 'source_name',
            doc => 'This column is deprectaed, do not use.',
            is_mutable => 0,
        },
    ],
    has_optional => [
        tgi_lims_status => { # TODO rename!
            is => 'Text',
            via => 'attributes',
            to => 'attribute_value',
            is_mutable => 1,
            is_many => 0,
            where => [ attribute_label => 'tgi_lims_status' ],
        },
        analysis_project_id => {
            is => 'Text',
            via => 'attributes',
            to => 'attribute_value',
            is_mutable => 1,
            is_many => 0,
            where => [ attribute_label => 'analysis_project_id' ],
        },
        original_est_fragment_size => {
            is => 'Integer',
            via => 'attributes',
            to => 'attribute_value',
            is_mutable => 1,
            is_many => 0,
            where => [ attribute_label => 'original_est_fragment_size' ],
        },
        original_est_fragment_size_max => {
            is => 'Integer',
            via => 'attributes',
            to => 'attribute_value',
            is_mutable => 1,
            is_many => 0,
            where => [ attribute_label => 'original_est_fragment_size_max' ],
        },
        original_est_fragment_size_min => {
            is => 'Integer',
            via => 'attributes',
            to => 'attribute_value',
            is_mutable => 1,
            is_many => 0,
            where => [ attribute_label => 'original_est_fragment_size_min' ],
        },
        original_est_fragment_std_dev => {
            is => 'Float',
            calculate_from => [ 'original_est_fragment_size_max', 'original_est_fragment_size_min', 'original_est_fragment_size' ],
            calculate => q(
                if ($original_est_fragment_size_max and
                    $original_est_fragment_size_min) {
                ($original_est_fragment_size_max - $original_est_fragment_size_min)/6; }
                else {
                    $original_est_fragment_size*.1;
                }),
        },
        read_orientation => {
            is => 'Text',
            via => 'attributes',
            to => 'attribute_value',
            is_mutable => 1,
            is_many => 0,
            where => [ attribute_label => 'read_orientation' ],
            valid_values => [ "forward_reverse", "reverse_forward" ],
        },
        run_name => {
            is => 'Text',
            len => 64,
        },
        subset_name => {
            is => 'Text',
            len => 64,
        },
        full_name => {
            calculate_from => [ 'run_name', 'subset_name' ],
            calculate => q("$run_name/$subset_name"),
        },
        full_path => {
            is => 'Text',
            via => 'attributes',
            to => 'attribute_value',
            is_mutable => 1,
            is_many => 0,
            where => [ attribute_label => 'full_path' ],
        },
        ignored => {
            is => 'Number',
            default_value => 0,
            via => 'attributes',
            to => 'attribute_value',
            is_mutable => 1,
            is_many => 0,
            where => [ attribute_label => 'ignored' ],
        },
        sample_source => {
            is => 'Genome::SampleSource',
            via => 'sample',
            to => 'source',
        },
        sample_source_id => {
            via => 'sample_source',
            to => 'id',
        },
        sample_source_name => {
            via => 'sample_source',
            to => 'name',
        },
        taxon => {
            is => 'Genome::Taxon',
            via => 'sample',
        },
        species_name => { via => 'taxon' },
    ],
    has_many_optional => [
        attributes => {
            is => 'Genome::InstrumentDataAttribute',
            reverse_as => 'instrument_data',
        },
        events => {
            is => 'Genome::Model::Event',
            reverse_as => 'instrument_data',
        },
        allocations => {
            is => 'Genome::Disk::Allocation',
            reverse_as => 'owner',
        },
    ],
    schema_name => 'GMSchema',
    data_source => 'Genome::DataSource::GMSchema',
    doc => 'Contains information common to all types of instrument data',
};

sub _preprocess_subclass_description {
    my ($class, $desc) = @_;

    for my $prop_name ( keys %{$desc->{has}} ) {
        my $prop_desc = $desc->{has}{$prop_name};
        next if not $prop_desc->{is_attribute};
        $prop_desc->{via} = 'attributes';
        $prop_desc->{where} = [ attribute_label => $prop_name ];
        $prop_desc->{to} = 'attribute_value';
        $prop_desc->{is_delegated} = 1;
        $prop_desc->{is_mutable} = 1;
    }

    return $desc;
}

sub create {
    my ($class) = @_;
    if ($class eq __PACKAGE__ or $class->__meta__->is_abstract) {
        # this class is abstract, and the super-class re-calls the constructor from the correct subclass
        return $class->SUPER::create(@_);
    }

    # This extra processing allows for someone to create instrument data with properties that aren't listed in any of the
    # class definitions. Instead of having UR catch these extras and die, they are captured here and later turned into
    # instrument data attributes.
    $class = shift;
    my ($bx, @extra);
    eval{ ($bx, @extra) = $class->define_boolexpr(@_); };
    return if not $bx;
    if ( @extra and @extra % 2 == 1 ) {
        $class->error_message("Odd number of attributes sent to create intrument data: ".Data::Dumper::Dumper(\@extra));
        return;
    }
    my %extra = @extra;

    my $self = $class->SUPER::create($bx);
    return if not $self;

    for my $key (keys %extra) {
        next unless length($key);
        next unless defined $extra{$key}; #skip undefined values
        my $attribute = Genome::InstrumentDataAttribute->create(
            attribute_label => $key,
            attribute_value => $extra{$key},
            instrument_data_id => $self->id,
        );
        unless ($attribute) {
            $self->error_message("Could not create attribute ".$key." => " . $extra{$key} . " for instrument data " . $self->id);
            $self->delete;
            return;
        }
    }

    return $self;
}

sub delete {
    my $self = shift;

    my ($expunge_status) = $self->_expunge_assignments;
    return unless $expunge_status;

    #finally, clean up the instrument data
    for my $attr ( $self->attributes ) {
        $attr->delete;
    }

    $self->_create_deallocate_observer;

    for my $attribute ($self->attributes) {
        $attribute->delete;
    }

    return $self->SUPER::delete;
}

my @_AR_SUBCLASSES_TO_EXPUNGE = (
    'Genome::InstrumentData::AlignmentResult::Merged',
    'Genome::InstrumentData::AlignmentResult',
    'Genome::InstrumentData::AlignmentResult::Tophat',
);

sub _expunge_assignments{
    my $self = shift;
    my $instrument_data_id = $self->id;
    my %affected_users;

    my @inputs = Genome::Model::Input->get(
        name => 'instrument_data',
        value_id => $instrument_data_id
    );
    my @models = map( $_->model, @inputs);

    for my $model (@models) {
        $model->remove_instrument_data($self);
        my $display_name = $self->__display_name__;
        push(@{$affected_users{$model->user_name}->{join(" ", $display_name, $self->id)}}, $model->id);
    }

    # There may be builds using this instrument data even though it had previously been unassigned from the model
    my @build_inputs = Genome::Model::Build::Input->get(
        name => 'instrument_data',
        value_id => $instrument_data_id
    );
    my @builds = map($_->build, @build_inputs);
    for my $build (@builds) {
        $build->abandon();
        push @models, $build->model;
    }

    for my $ar_subclass_name (@_AR_SUBCLASSES_TO_EXPUNGE) {
        my @alignment_results = $ar_subclass_name->get(
            instrument_data_id => $self->id);

        for my $alignment_result (@alignment_results) {
            for my $bamqc_result (Genome::InstrumentData::AlignmentResult::Merged::BamQc->get(alignment_result_id => $alignment_result->id)) {
                unless ($bamqc_result->delete) {
                    die $self->error_message(sprintf(
                            "Failed to remove BamQc result (%s) for alignment result (%s)",
                            $alignment_result->id, $bamqc_result->id));
                }
            }

            unless ($alignment_result->delete) {
                die $self->error_message("Could not remove instrument data " . $self->__display_name__ .
                    " because alignment result " . $alignment_result->__display_name__ .
                    " that uses this instrument data could not be deleted!");
            }
        }

    }

    return 1, %affected_users;
}

sub _create_deallocate_observer {
    my $self = shift;
    my @allocations = $self->allocations;
    return 1 unless @allocations;
    my $deallocator;
    $deallocator = sub {
        for my $allocation (@allocations) {
            $allocation->delete;
        }
        UR::Context->cancel_change_subscription(
            'commit', $deallocator
        );
    };
    UR::Context->create_subscription(
        method => 'commit',
        callback => $deallocator
    );
    return 1;
}

sub calculate_alignment_estimated_kb_usage {
    my $self = shift;
    Carp::confess "calculate_alignment_estimated_kb_usage not overridden in instrument data subclass " . $self->class;
}

sub sample_type {
    my $self = shift;
    my $sample_extraction_type = $self->sample->extraction_type;
    return unless defined $sample_extraction_type;

    if ($sample_extraction_type eq 'genomic_dna' or $sample_extraction_type eq 'pooled dna') {
        return 'dna';
    }
    elsif ($sample_extraction_type eq 'rna') {
        return 'rna';
    }
    return;
}

sub create_mock {
    my $class = shift;
    return $class->SUPER::create_mock(subclass_name => 'Genome::InstrumentData', @_);
}

sub run_identifier  {
    die "run_identifier not defined in instrument data subclass.  please define this. this method should " .
         "provide a unique identifier for the experiment/run (eg flow_cell_id, ptp barcode, etc).";
}

sub dump_fastqs_from_bam {
    my $self = shift;
    my %p = @_;

    unless ($self->can('bam_path')) {
        my $error = 'cannot call bam path';
        $self->error_message($error);
        die $error;
    }
    unless ($self->bam_path) {
        my $error = 'Attempted to dump a bam but the path is not defined.';
        $self->error_message($error);
        die $error;
    }
    unless (-e $self->bam_path) {
        my $error = 'Attempted to dump a bam but the path does not exist:' . $self->bam_path;
        $self->error_message($error);
        die $error;
    }

    my $directory = delete $p{directory};
    if ( not $directory ) {
        $directory = Genome::Sys->base_temp_directory.'/'.'unpacked_bam_'.$self->id;
        Genome::Sys->create_directory($directory) if not -d $directory;
    }

    my $subset = (defined $self->subset_name ? $self->subset_name : 0);

    my %read_group_params;

    if (defined $p{read_group_id}) {
        $read_group_params{read_group_id} = delete $p{read_group_id};
        $self->status_message("Using read group id " . $read_group_params{read_group_id});
    }

    my $fwd_file = sprintf("%s/s_%s_1_sequence.txt", $directory, $subset);
    my $rev_file = sprintf("%s/s_%s_2_sequence.txt", $directory, $subset);
    my $fragment_file = sprintf("%s/s_%s_sequence.txt", $directory, $subset);
    my $cmd = Genome::Model::Tools::Picard::SamToFastq->create(input=>$self->bam_path, fastq=>$fwd_file, fastq2=>$rev_file, fragment_fastq=>$fragment_file, no_orphans=>1, %read_group_params);
    if ( not $cmd ) {
        die $self->error_message('Failed to create gmt picard sam-to-fastq');
    }
    $cmd->dump_status_messages(1);
    unless ($cmd->execute()) {
        die $cmd->error_message;
    }

    if ((-s $fwd_file && !-s $rev_file) ||
        (!-s $fwd_file && -s $rev_file)) {
        $self->error_message("Fwd & Rev files are lopsided; one has content and the other doesn't. Can't proceed");
        die $self->error_message;
    }

    my @files;
    if (-s $fwd_file && -s $rev_file) {
        push @files, ($fwd_file, $rev_file);
    }
    if (-s $fragment_file && !$p{discard_fragments}) {
        push @files, $fragment_file;
    }

    return @files;
}

sub lane_qc_models {
    my $self = shift;

    # Find the Lane QC models that use this instrument data
    my $instrument_data_id = $self->id;
    my @inputs = Genome::Model::Input->get(value_id => $instrument_data_id);
    my @lane_qc_models = grep { $_->is_lane_qc }
                         grep { $_->class eq 'Genome::Model::ReferenceAlignment' }
                         map  { $_->model } @inputs;

    # Find the Lane QC models that used the default genotype_microarray input
    my $sample = $self->sample;
    my @gm_model_ids = map { $_->id } ($sample->default_genotype_models);
    @lane_qc_models = grep { $_->inputs(name => 'genotype_microarray', value_id => \@gm_model_ids) } @lane_qc_models;

    return @lane_qc_models;;
}

sub lane_qc_build {
    my $self = shift;
    my @qc_models = $self->lane_qc_models;
    return unless @qc_models;
    # find the newest succeeded build
    my @builds = sort { $b->date_scheduled cmp $a->date_scheduled }
                 map { $_->succeeded_builds }
                 @qc_models;
    return unless @builds;
    return $builds[0];
}

sub lane_qc_dir {
    my $self = shift;
    my $qc_build = $self->lane_qc_build;
    return unless ($qc_build);
    my $qc_dir = $qc_build->qc_directory;
    if (-d $qc_dir) {
        return $qc_dir;
    }
    else {
        return;
    }
}


sub __display_name__ {
    my $self = shift;
    my $subset_name = $self->subset_name || 'unknown-subset';
    my $run_name = $self->run_name || 'unknown-run-name';
    return join '.', $run_name, $subset_name;
}

1;

__END__
s:
_sls_fwd_read_length
_sls_read_length
_sls_rev_read_length
adaptor_path	attributes
analysis_software_version	attributes
archive_path	attributes
bam_path	attributes
clusters	attributes
cycles
fastqc_path	attributes
filt_error_rate_avg
flow_cell
flow_cell_id	attributes
fwd_clusters	attributes
fwd_filt_aligned_clusters_pct	attributes
fwd_filt_error_rate_avg
fwd_kilobases_read	attributes
fwd_read_length	attributes
fwd_run_type	attributes
fwd_seq_id	attributes
gc_bias_path	attributes
gerald_directory	attributes
index_sequence	attributes
is_external	attributes
is_paired_end
lane	attributes
median_insert_size
old_filt_error_rate_avg	attributes
old_fwd_filt_error_rate_avg	attributes
old_median_insert_size	attributes
old_rev_filt_error_rate_avg	attributes
old_sd_above_insert_size	attributes
old_sd_below_insert_size	attributes
project
project_name	attributes
read_1_pct_mismatch
read_2_pct_mismatch
read_count
read_length	attributes
rev_clusters	attributes
rev_filt_aligned_clusters_pct	attributes
rev_filt_error_rate_avg
rev_kilobases_read	attributes
rev_read_length	attributes
rev_run_type	attributes
rev_seq_id	attributes
run_type	attributes
sd_above_insert_size
sd_below_insert_size
sd_insert_size
sequencing_platform
short_name
subclass_name
target_region_set_name	attributes

i:
barcode	attributes
base_count	attributes
blacklisted_segments	attributes
description	attributes
fragment_count	attributes
full_name
fwd_read_length	attributes
genotype_file	attributes
import_date	attributes
import_format	attributes
import_source_name	attributes
is_paired_end	attributes
median_insert_size	attributes
original_data_path	attributes
read_count	attributes
read_length	attributes
reference_sequence_build
reference_sequence_build_id	attributes
rev_read_length	attributes
sd_above_insert_size	attributes
source	sample
source_id	source
source_name	source
sra_accession	attributes
sra_sample_id	attributes
subclass_name
target_region_set_name	attributes
user_name	attributes

f:
_fasta_file
_qual_file
beads_loaded	attributes
copies_per_bead	attributes
fc_id	attributes
incoming_dna_name	attributes
index_sequence	attributes
is_paired_end	attributes
key_pass_wells	attributes
predicted_recovery_beads	attributes
read_count	attributes
region_id	attributes
region_number	attributes
research_project	attributes
sample_set	attributes
sequencing_platform
sff_file	attributes
ss_id	attributes
subclass_name
supernatant_beads	attributes
total_bases_read	attributes
total_key_pass	attributes
total_raw_wells	attributes
total_reads	attributes

g:
disk_allocation
read_count
research_project	attributes
sequencing_platform
subclass_name

