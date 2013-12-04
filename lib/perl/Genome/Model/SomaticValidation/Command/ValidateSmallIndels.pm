package Genome::Model::SomaticValidation::Command::ValidateSmallIndels;

use strict;
use warnings;

use Genome;

use File::Spec qw();
use File::Basename;
use Workflow::Simple;

class Genome::Model::SomaticValidation::Command::ValidateSmallIndels {
    is => 'Genome::Command::Base',
    has => [
        varscan_version => {
            is => 'Text',
            doc => "Varscan version to use when running varscan validation" ,
        },
        samtools_use_baq => {
            is => 'Boolean',
            doc => 'When doing pileup/mpileup, should we enable baq (-B) option',
            default_value => 0,
        },
    ],
    has_optional => [
        build => {
            is => 'Genome::Model::Build::SomaticValidation',
            id_by => 'build_id',
        },
        build_id => {
            is => 'Text',
            is_input => 1,
            doc => 'build id of SomaticValidation model. Provide this OR the other options.',
        },
        small_indel_output_bed => {
            is => 'Text',
            doc => "File of small indels to be realigned.",
        },
        tumor_bam   => {
            is => 'Text',
            doc => "Tumor Bam File (Validation Bam)",
            is_input => 1,
        },
        normal_bam  => {
            is => 'Text',
            doc => "Normal Bam File (Validation Bam)",
            is_input => 1,
        },
        reference_fasta => {
            is => 'Text',
            doc => "Reference Fasta" ,
            is_optional => 1,
        },
        varscan_indel_output => {
            is => 'Text',
            doc => "gmt varscan validation output run on realigned bams" ,
        },
        varscan_snp_output  => {
            is => 'Text',
            doc => "gmt varscan validation output run on realigned bams" ,
        },
        final_output_file   => {
            is => 'Text',
            doc => "gmt varscan process-validation-indels final output file labeling indels as Somatic or otherwise" ,
        },
        realigned_bam_file_directory => {
            is => 'Text',
            doc => "Where to dump the realigned bam file",
        },
        output_dir => {
            is => 'Text',
            doc => "Base output directory if the build is not set",
            is_input => 1,
        },
        # FIXME fix up these three params
        varscan_params => {
            calculate_from => [qw/ normal_purity min_var_frequency /],
            calculate => q| '--validation 1 --somatic-p-value 1.0e-02 --p-value 0.10 --min-coverage 8 --min-var-freq '.$min_var_frequency.' --normal-purity '.$normal_purity |,
        },
        normal_purity => {
            is => 'Float',
            doc => "Normal purity param to pass to varscan",
            default => 1
        },
        min_var_frequency => {
            is => 'Float',
            doc => "Minimum variant frequency to pass to varscan",
            default => 0.08
        },
    ],
    has_param => [
        lsf_queue => {
            default => $ENV{GENOME_LSF_QUEUE_BUILD_WORKER_ALT},
        },
    ],
};

sub execute {
    my $self = shift;

    return 1 if $self->build and not $self->build->normal_sample;

    $self->_resolve_inputs;

    $self->_create_output_directory();

    $self->_run_workflow;

    return 1;
}

sub _run_workflow {
    my $self = shift;

    # Define a workflow from the static XML at the bottom of this module
    my $workflow = Workflow::Operation->create_from_xml(\*DATA);
    # Validate the workflow
    my @errors = $workflow->validate;
    if (@errors) {
        $self->error_message(@errors);
        die "Errors validating workflow\n";
    }

    # Collect and set input parameters
    # TODO stop hardcoding stuff
    my %input;

    # Params for gatk
    $input{gatk_memory} = '16';
    $input{gatk_version} = 5777;
    $input{index_bam} = 1,
    $input{target_intervals_are_sorted} = 0;
    $input{small_indel_list} = $self->small_indel_output_bed;
    $input{normal_bam} = $self->normal_bam;
    $input{tumor_bam} = $self->tumor_bam;
    $input{realigned_tumor_bam} = $self->_realigned_tumor_bam_file;
    $input{realigned_normal_bam} = $self->_realigned_normal_bam_file;
    $input{reference} = $self->reference_fasta;

    # Params for varscan
    $input{samtools_use_baq} = $self->samtools_use_baq;
    $input{varscan_version} = $self->varscan_version;
    $input{varscan_indel_output} = $self->varscan_indel_output;
    $input{varscan_snp_output}= $self->varscan_snp_output;
    $input{varscan_params} = $self->varscan_params;
    my $bed = $self->small_indel_output_bed;
    ($input{small_indel_list_nobed} = $bed) =~ s/\.padded1bp\.bed$/\.annotation_format/;
    $input{final_output_file} = $self->final_output_file;

    my $log_dir = $self->output_dir;
    if(Workflow::Model->parent_workflow_log_dir) {
        $log_dir = Workflow::Model->parent_workflow_log_dir;
    }

    $workflow->log_dir($log_dir);

    # Launch workflow
    $self->status_message("Launching workflow now.");
    my $result = Workflow::Simple::run_workflow_lsf( $workflow, %input);

    # Collect and analyze results
    unless($result){
        die $self->error_message("Workflow did not return correctly.");
    }

    return 1;
}

sub _realigned_normal_bam_file {
    my $self = shift;
    my $realigned_normal_bam_file = basename($self->normal_bam,qr{\.bam});
    return $self->realigned_bam_file_directory . "/$realigned_normal_bam_file.realigned.bam";
}

sub _realigned_tumor_bam_file {
    my $self = shift;
    my $realigned_tumor_bam_file = basename($self->tumor_bam,qr{\.bam});
    return $self->realigned_bam_file_directory . "/$realigned_tumor_bam_file.realigned.bam";
}

sub _resolve_output_directory {
    my $self = shift;
    my $output_dir;
    if ($self->build) {
        $output_dir = File::Spec->join($self->build->data_directory, 'validation/small_indel');
    }
    return $output_dir;
}

sub _create_output_directory {
    my $self = shift;
    my $output_directory = $self->_resolve_output_directory();
    if ($output_directory) {
        Genome::Sys->create_directory($output_directory);
    }
    return $output_directory;
}

# If a build is provided, use that to populate other inputs. Otherwise make sure the other inputs are set. Primarily this is to make this testable.
sub _resolve_inputs {
    my $self = shift;

    # Make sure that if output paths arent set that the somatic variation build is, and set good defaults
    if ($self->build) {
        if ($self->final_output_file || $self->realigned_bam_file_directory || $self->small_indel_output_bed
            || $self->varscan_indel_output || $self->varscan_snp_output) {
            die $self->error_message("If a build is provided, you should not provide other params");
        }

        my $build = $self->build;
        my $model= $build->model;
        my $base_dir = $self->_resolve_output_directory();
        $self->final_output_file("$base_dir/final_output");
        $self->realigned_bam_file_directory("$base_dir/realigned_bams");
        $self->small_indel_output_bed("$base_dir/small_indels.bed");
        $self->varscan_indel_output("$base_dir/varscan_indels");
        $self->varscan_snp_output("$base_dir/varscan_snps");
        $self->tumor_bam($build->tumor_bam);
        $self->normal_bam($build->normal_bam);
        my $ref_seq_build_id = $model->reference_sequence_build->build_id;
        my $ref_seq_build = Genome::Model::Build->get($ref_seq_build_id);
        my $reference = $ref_seq_build->full_consensus_path('fa');
        $self->reference_fasta($reference);
    } else {
        my @required_properties = qw(final_output_file realigned_bam_file_directory small_indel_output_bed varscan_indel_output varscan_snp_output tumor_bam normal_bam reference_fasta);
        my $fail = 0;
        for my $property (@required_properties) {
            unless (defined $self->$property) {
                $fail = 1;
                $self->error_message("$property is not set and must be if somatic_validation_build is not set");
            }
        }
        die $self->error_message("All of the above properties must be set unless somatic_validation_build is set.") if $fail;
    }

    return 1;
}

1;

__DATA__
<?xml version='1.0' standalone='yes'?>

<workflow name="Validate Small Indels Subworkflow">

  <link fromOperation="input connector" fromProperty="gatk_memory" toOperation="Realign Indels Tumor" toProperty="max_memory" />
  <link fromOperation="input connector" fromProperty="gatk_version" toOperation="Realign Indels Tumor" toProperty="version" />
  <link fromOperation="input connector" fromProperty="small_indel_list" toOperation="Realign Indels Tumor" toProperty="target_intervals" />
  <link fromOperation="input connector" fromProperty="realigned_tumor_bam" toOperation="Realign Indels Tumor" toProperty="output_realigned_bam" />
  <link fromOperation="input connector" fromProperty="tumor_bam" toOperation="Realign Indels Tumor" toProperty="input_bam" />
  <link fromOperation="input connector" fromProperty="reference" toOperation="Realign Indels Tumor" toProperty="reference_fasta" />
  <link fromOperation="input connector" fromProperty="target_intervals_are_sorted" toOperation="Realign Indels Tumor" toProperty="target_intervals_are_sorted" />
  <link fromOperation="input connector" fromProperty="index_bam" toOperation="Realign Indels Tumor" toProperty="index_bam" />

  <link fromOperation="input connector" fromProperty="gatk_memory" toOperation="Realign Indels Normal" toProperty="max_memory" />
  <link fromOperation="input connector" fromProperty="gatk_version" toOperation="Realign Indels Normal" toProperty="version" />
  <link fromOperation="input connector" fromProperty="small_indel_list" toOperation="Realign Indels Normal" toProperty="target_intervals" />
  <link fromOperation="input connector" fromProperty="realigned_normal_bam" toOperation="Realign Indels Normal" toProperty="output_realigned_bam" />
  <link fromOperation="input connector" fromProperty="normal_bam" toOperation="Realign Indels Normal" toProperty="input_bam" />
  <link fromOperation="input connector" fromProperty="reference" toOperation="Realign Indels Normal" toProperty="reference_fasta" />
  <link fromOperation="input connector" fromProperty="target_intervals_are_sorted" toOperation="Realign Indels Normal" toProperty="target_intervals_are_sorted" />
  <link fromOperation="input connector" fromProperty="index_bam" toOperation="Realign Indels Normal" toProperty="index_bam" />

  <link fromOperation="input connector" fromProperty="samtools_use_baq" toOperation="Varscan Validation" toProperty="samtools_use_baq" />
  <link fromOperation="input connector" fromProperty="reference" toOperation="Varscan Validation" toProperty="reference" />
  <link fromOperation="input connector" fromProperty="varscan_version" toOperation="Varscan Validation" toProperty="version" />
  <link fromOperation="input connector" fromProperty="varscan_indel_output" toOperation="Varscan Validation" toProperty="output_indel" />
  <link fromOperation="input connector" fromProperty="varscan_snp_output" toOperation="Varscan Validation" toProperty="output_snp" />
  <link fromOperation="input connector" fromProperty="varscan_params" toOperation="Varscan Validation" toProperty="varscan_params" />
  <link fromOperation="Realign Indels Tumor" fromProperty="output_realigned_bam" toOperation="Varscan Validation" toProperty="tumor_bam" />
  <link fromOperation="Realign Indels Normal" fromProperty="output_realigned_bam" toOperation="Varscan Validation" toProperty="normal_bam" />

  <link fromOperation="input connector" fromProperty="final_output_file" toOperation="Process Validation Indels" toProperty="output_file" />
  <link fromOperation="input connector" fromProperty="small_indel_list_nobed" toOperation="Process Validation Indels" toProperty="variants_file" />
  <link fromOperation="Varscan Validation" fromProperty="output_indel" toOperation="Process Validation Indels" toProperty="validation_indel_file" />
  <link fromOperation="Varscan Validation" fromProperty="output_snp" toOperation="Process Validation Indels" toProperty="validation_snp_file" />

  <link fromOperation="Process Validation Indels" fromProperty="output_file" toOperation="output connector" toProperty="output_file" />

  <operation name="Realign Indels Tumor">
    <operationtype commandClass="Genome::Model::Tools::Gatk::IndelRealigner" typeClass="Workflow::OperationType::Command" />
  </operation>
  <operation name="Realign Indels Normal">
    <operationtype commandClass="Genome::Model::Tools::Gatk::IndelRealigner" typeClass="Workflow::OperationType::Command" />
  </operation>
  <operation name="Varscan Validation">
    <operationtype commandClass="Genome::Model::Tools::Varscan::Validation" typeClass="Workflow::OperationType::Command" />
  </operation>
  <operation name="Process Validation Indels">
    <operationtype commandClass="Genome::Model::Tools::Varscan::ProcessValidationIndels" typeClass="Workflow::OperationType::Command" />
  </operation>

  <operationtype typeClass="Workflow::OperationType::Model">
    <inputproperty>samtools_use_baq</inputproperty>
    <inputproperty>normal_bam</inputproperty>
    <inputproperty>gatk_memory</inputproperty>
    <inputproperty>gatk_version</inputproperty>
    <inputproperty>index_bam</inputproperty>
    <inputproperty>target_intervals_are_sorted</inputproperty>
    <inputproperty>small_indel_list</inputproperty>
    <inputproperty>normal_bam</inputproperty>
    <inputproperty>tumor_bam</inputproperty>
    <inputproperty>realigned_tumor_bam</inputproperty>
    <inputproperty>realigned_normal_bam</inputproperty>
    <inputproperty>reference</inputproperty>
    <inputproperty>varscan_version</inputproperty>
    <inputproperty>varscan_indel_output</inputproperty>
    <inputproperty>varscan_snp_output</inputproperty>
    <inputproperty>varscan_params</inputproperty>
    <inputproperty>small_indel_list_nobed</inputproperty>
    <inputproperty>final_output_file</inputproperty>

    <outputproperty>output_file</outputproperty>
  </operationtype>

</workflow>
