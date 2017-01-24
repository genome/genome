package Genome::InstrumentData;

use strict;
use warnings;

use Genome;

use List::MoreUtils qw(uniq);

class Genome::InstrumentData {
    roles => ['Genome::Role::ObjectWithAllocations', 'Genome::Role::Notable'],
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
            to => 'individual',
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
        analysis_project_bridges => {
            is => 'Genome::Config::AnalysisProject::InstrumentDataBridge',
            reverse_as => 'instrument_data',
            is_many => 1,
        },
        analysis_projects => {
            is => 'Genome::Config::AnalysisProject',
            via => 'analysis_project_bridges',
            to => 'analysis_project',
            is_many => 1,
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
        $class->error_message("Odd number of attributes sent to create instrument data: ".Data::Dumper::Dumper(\@extra));
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

    for my $attribute ($self->attributes) {
        $attribute->delete;
    }

    for my $bridge ($self->analysis_project_bridges) {
        $bridge->delete;
    }

    return $self->SUPER::delete;
}

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
        push(@{$affected_users{$model->created_by}->{join(" ", $display_name, $self->id)}}, $model->id);
        push(@{$affected_users{$model->run_as}->{join(" ", $display_name, $self->id)}}, $model->id);
    }

    # There may be builds using this instrument data even though it had previously been unassigned from the model
    my @build_inputs = Genome::Model::Build::Input->get(
        name => 'instrument_data',
        value_id => $instrument_data_id
    );
    my @builds = map($_->build, @build_inputs);
    for my $build (@builds) {
        $build->abandon(undef, 'Expunging instrument data ' . $instrument_data_id);
        push @models, $build->model;
    }

    Genome::SoftwareResult->expunge_results_containing_object($self, 'Expunging instrument data ' . $instrument_data_id);

    return 1, %affected_users;
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

sub alignment_results_from_analysis_projects {
    my $self = shift;

    my @aps = $self->analysis_projects;
    return unless scalar(@aps);

    my @inputs = Genome::Model::Input->get(value_id => $self->id, 'model.analysis_projects.id' => [map $_->id, @aps]);
    return unless scalar(@inputs);

    my @alignment_results;
    for my $input (@inputs) {
        my $build = $input->model->last_complete_build;
        next unless $build;
        next unless Genome::Model::Build::Input->get(value_id => $self->id, build_id => $build->id);

        push @alignment_results,
            grep { $_->instrument_data_id eq $self->id }
            grep { $_->isa('Genome::InstrumentData::AlignmentResult') }
            $build->all_results;
    }

    return uniq(@alignment_results);
}

sub projects {
    my $self = shift;

    my @parts = Genome::ProjectPart->get(
        entity_id => $self->id,
        entity_class_name=> 'Genome::InstrumentData',
        label => 'instrument_data',
    );
    my @projects = Genome::Project->get([map $_->project_id, @parts]);
    return @projects;
}

sub __display_name__ {
    my $self = shift;
    my $subset_name = $self->subset_name || 'unknown-subset';
    my $run_name = $self->run_name || 'unknown-run-name';
    return join '.', $run_name, $subset_name;
}

1;
