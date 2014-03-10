package Genome::Model::Build::SomaticValidation::IdentifyDnpResult;

class Genome::Model::Build::SomaticValidation::IdentifyDnpResult {
    is => 'Genome::SoftwareResult::Stageable',
    has_param => [
        proportion => {
            is => 'Number',
            doc => 'Proportion of reads supporting the DNP required for the site to be considered a DNP.',
        },
    ],
    has_input => [
        dv2_result_id => {
            is => 'Text',
        },
        tumor_aligment_result_id => {
            is => 'Text',
        },
    ],
    has_optional => [
        dv2_result => {
            is => 'Genome::Model::Tools::DetectVariants2::Result::Base',
            id_by => 'dv2_result_id',
        },
        tumor_alignment_result => {
            is => 'Genome::InstrumentData::AlignedBamResult',
            id_by => 'tumor_aligment_result_id',
        },
    ],
};

sub create {
    my $class = shift;

    my $self = $class->SUPER::create(@_);
    return unless $self;

    my $rv = eval {
        $self->_prepare_staging_directory;
        $self->_generate_data;
        $self->_prepare_output_directory;
        $self->_promote_data;
        $self->_reallocate_disk_allocation;
        return 1;
    };
    my $error = $@;

    if ($error) {
        die $self->error_message($error);
    }
    elsif ($rv ne 1) {
        die $self->error_message("Unexpected return value: $rv");
    }

    $self->debug_message('All processes completed.');

    return $self;
}

sub _snvs_hq_bed {
    my $self = shift;
    my $version = 2;
    my $snvs_hq_bed = join('/', $self->dv2_result->output_dir, 'snvs.hq.bed');
    unless (defined $snvs_hq_bed) {
        die $self->error_message("Could not resolve snvs.hq.bed from input build.");
    }
    return $snvs_hq_bed;
}

sub _tumor_bam_file {
    my $self = shift;
    my $bam_file = $self->tumor_alignment_result->bam_path;
    unless (defined $bam_file) {
        die $self->error_message("'whole_rmdup_bam_file' not defined for build's tumor_reference_alignment.");
    }
    return $bam_file;
}

sub _generate_data {
    my $self = shift;

    my $annotation_file = Genome::Sys->create_temp_file_path;
    my $snvs_hq_bed = $self->_snvs_hq_bed;
    my $bed_to_annotation = Genome::Model::Tools::Bed::Convert::BedToAnnotation->create(
        snv_file => $snvs_hq_bed,
        output => $annotation_file,
    );
    unless($bed_to_annotation->execute) {
        die $self->error_message('Failed to convert snvs.hq.bed to annotation format.');
    }

    my $bam_file = $self->_tumor_bam_file;
    my $output_file = join('/', $self->temp_staging_directory, 'identify_dnp_output.tsv');

    #FIXME make sure annotation_file is sorted (chr and coord)
    my $identify_dnp = Genome::Model::Tools::Somatic::IdentifyDnp->create(
        annotation_input_file => $annotation_file,
        proportion => $self->proportion,
        bam_file => $bam_file,
        output_file => $output_file,
    );
    unless($identify_dnp->execute) {
        die $self->error_message('Failed to complete "Identify DNP" command.');
    }

    return 1;
}

sub required_rusage { '' };
sub _needs_symlinks_followed_when_syncing { 0 };
sub _working_dir_prefix { 'somatic-validation-result' };
sub resolve_allocation_disk_group_name { $ENV{GENOME_DISK_GROUP_MODELS} };

sub resolve_allocation_subdirectory {
    my $self = shift;
    my $staged_basename = File::Basename::basename($self->temp_staging_directory);
    return join('/', 'build_merged_alignments', $self->id, 'identify-dnp-' . $staged_basename);
};

sub estimated_kb_usage {
    my $self = shift;
    my $snvs_hq_bed = $self->_snvs_hq_bed;
    my $bytes = -s $snvs_hq_bed;
    return (2 * $bytes / 1000);
}

1;
