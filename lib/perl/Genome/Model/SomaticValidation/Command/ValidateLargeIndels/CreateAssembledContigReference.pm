package Genome::Model::SomaticValidation::Command::ValidateLargeIndels::CreateAssembledContigReference;

use strict;
use warnings;
use Genome;

use Data::UUID qw();
use File::Spec qw();

class Genome::Model::SomaticValidation::Command::ValidateLargeIndels::CreateAssembledContigReference {
    is => 'Command::V2',
    has => [
        build => {
            is => 'Genome::Model::Build::SomaticValidation',
            id_by => 'build_id',
        },
        build_id => {
            is => 'Text',
            is_input => 1,
            is_output => 1,
            doc => 'build id of SomaticValidation model',
        },
        _long_indel_bed_file => {
            is => 'Text',
            is_input => 1,
            is_optional => 1,
            doc => 'Path to bed file that contains the large indels. Resolved via the build.',
        },
        reference_transcripts => {
            is => 'Text',
            default => "NCBI-human.ensembl/67_37l_v2", #TODO this should be a param from the somatic validation processing profile
            #example_values => ["NCBI-human.ensembl/67_37l_v2"], #TODO this should be a param from the somatic validation processing profile
            doc => 'The set of reference_transcripts to use, which get passed to the annotator',
        },
    ],
    has_output => [
        reference_build_id => {
            is => 'Text',
            doc => 'ID of newly created reference sequence',
            is_optional => 1,
        },
        skip => {
            is => 'Boolean',
            doc => 'indicate whether large indel validation should take place',
            is_optional => 1,
        },
        output_dir => {
            is => 'Text',
            doc => 'Location of output directory',
            is_optional => 1,
        },
        contigs_fasta => {
            is => 'Text',
            doc => 'Location of the new contigs',
            is_optional => 1,
        },
        tier_files => {
            is => 'Text',
            doc => 'Location of the tier files',
            is_optional => 1,
        },
    ],
    has_param => [
        lsf_queue => {
            default => Genome::Config::get('lsf_queue_build_worker_alt'),
        },
    ],
};

sub sub_command_category {'pipeline steps'}

sub shortcut {
    my $self = shift;

    return $self->skip_validation;
}

sub execute {
    my $self = shift;

    return 1 if $self->skip_validation;

    my $output_directory = $self->_create_output_directory;

    $self->reference_build_id(1);
    $self->contigs_fasta(1);
    $self->tier_files(1);

    my $long_indel_bed_file = $self->_resolve_long_indel_bed_file;
    $self->_long_indel_bed_file($long_indel_bed_file);

    my $sample_id = Data::UUID->new->create_str();
    #TODO the instructions say to "be sure to save the STDOUT from this tool

    my $ref_seq_build = $self->build->reference_sequence_build;
    my $ref_seq_fasta = $ref_seq_build->full_consensus_path('fa');
    my $annotator_version = $self->build->processing_profile->transcript_variant_annotator_version;
    my $tumor_bam = $self->build->tumor_bam;
    my $normal_bam = $self->build->normal_bam;

    my $cmd = Genome::Model::Tools::Validation::LongIndelsGenerateMergedAssemblies->create(
        long_indel_bed_file => $self->_long_indel_bed_file,
        output_dir => $output_directory,
        transcript_variant_annotator_version => $annotator_version,
        reference_transcripts => $self->reference_transcripts,
        tumor_bam => $tumor_bam,
        normal_bam => $normal_bam,
        reference_fasta => $ref_seq_fasta,
    );

    unless ($cmd->execute) {
        die $self->error_message("Failed to generate merged assemblies");
    }

    my $contigs_file = $cmd->contigs_fasta;
    unless (-s $contigs_file) {
        $self->warning_message("Failed to get valid assembly contig fasta. Skipping remainder of long indel validation.");
        $self->skip(1);
        return 1;
    }

    $self->contigs_fasta($contigs_file);

    my $annotation_build = Genome::Model::Build::ImportedAnnotation->get(name => $self->reference_transcripts);
    my $tier_file_location = $annotation_build->tiering_bed_files_by_version(3);
    $self->tier_files($tier_file_location);

    my %import_params = (
        species_name => 'human',
        use_default_sequence_uri => '1',
        derived_from => $ref_seq_build,
        append_to => $ref_seq_build,
        version => '500bp_assembled_contigs',
        fasta_file => $contigs_file,
        prefix => $sample_id,
        is_rederivable => 1,
    );
    if ($self->build->model->analysis_project) {
        $import_params{analysis_project} = $self->build->model->analysis_project;
    }

    my $new_ref_cmd = Genome::Model::Command::Define::ImportedReferenceSequence->create(%import_params);
    unless ($new_ref_cmd->execute) {
        $self->error_message('Failed to execute the definition of the new reference sequence with added contigs.');
        return;
    }

    my $new_ref_build_id = $new_ref_cmd->result_build_id;
    my $new_ref_build = Genome::Model::Build->get($new_ref_build_id);
    unless ($new_ref_build->status eq 'Succeeded') {
        die $self->error_message("New reference build not successful.");
    }

    $self->reference_build_id($new_ref_build_id);

    return 1;
}

sub _resolve_long_indel_bed_file {
    my $self = shift;
    my $long_indel_bed_file = $self->build->data_directory . "/validation/small_indel/large_indels.bed";
    unless (-s $long_indel_bed_file) {
        $self->warning_message("Long indel bed file $long_indel_bed_file does not exist or has no size");
        return;
    }
    return $long_indel_bed_file;
}

sub _resolve_output_directory {
    my $self = shift;
    my $output_dir = File::Spec->join($self->build->data_directory, '/validation/large_indel');
    $self->output_dir($output_dir);
    return $output_dir;
}

sub _create_output_directory {
    my $self = shift;

    my $output_directory = $self->_resolve_output_directory();
    Genome::Sys->create_directory($output_directory);

    return $output_directory;
}

sub skip_validation {
    my $self = shift;

    my $should_skip = 0;
    unless($self->build->run_indel_validation) {
        $self->debug_message('Build indicates indel validation should be skipped.');
        $should_skip = 1;
    }

    my $long_indel_bed_file = $self->_resolve_long_indel_bed_file;
    my $skip_msg = 'Skipping long indel validation';

    unless ($long_indel_bed_file) {
        $self->warning_message("No long indel bed file exists with size. $skip_msg");
        $should_skip = 1;
    }
    unless ($self->build->normal_sample) {
        $self->warning_message("Somatic validation of a single bam.  $skip_msg");
        $should_skip = 1;
    }

    if ($should_skip) {
        $self->skip(1);
        $self->_resolve_output_directory;

        $self->reference_build_id(1);
        $self->contigs_fasta(1);
        $self->tier_files(1);
    } else {
        $self->skip(0);
    }

    return $self->skip;
}

1;

