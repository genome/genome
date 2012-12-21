package Genome::Model::SomaticValidation::Command::ValidateSmallIndels;

use strict;
use warnings;

use Genome;

use File::Spec qw();
use File::Basename;

class Genome::Model::SomaticValidation::Command::ValidateSmallIndels {
    is => 'Genome::Command::Base',
    has_optional => [
        build => {
            is => 'Genome::Model::Build::SomaticValidation',
            id_by => 'build_id',
        },
        build_id => {
            is => 'Integer',
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
        },
        normal_bam  => {
            is => 'Text',
            doc => "Normal Bam File (Validation Bam)",
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
            default => 'apipe',
        },
    ],
};

sub execute {
    my $self = shift;

    $self->_resolve_inputs;

    $self->_create_output_directory();

    $self->_run_gatk;

    $self->_run_varscan;

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

sub _run_gatk {
    my $self = shift;

    my $small_indel_list = $self->small_indel_output_bed;
    my $normal_bam = $self->normal_bam;
    my $tumor_bam = $self->tumor_bam;
    my $realigned_tumor_bam_file = $self->_realigned_tumor_bam_file;
    my $realigned_normal_bam_file = $self->_realigned_normal_bam_file;
    my $reference = $self->reference_fasta;

    my $gatk_tumor_cmd = Genome::Model::Tools::Gatk::RealignIndels->create(
        max_memory => "16g",
        version => 5777,
        target_intervals => $small_indel_list,
        output_realigned_bam => $realigned_tumor_bam_file,
        input_bam => $tumor_bam,
        reference_fasta => $reference,
        target_intervals_are_sorted => 0,
    );

    my $gatk_normal_cmd = Genome::Model::Tools::Gatk::RealignIndels->create(
        max_memory => "16g",
        version => 5777,
        target_intervals => $small_indel_list,
        output_realigned_bam => $realigned_normal_bam_file,
        input_bam => $normal_bam,
        reference_fasta => $reference,
        target_intervals_are_sorted => 0,
    );

    unless ($gatk_tumor_cmd->execute) {
        die $self->error_message("Failed to run gatk IndelRealigner on tumor");
    }

    unless ($gatk_normal_cmd->execute) {
        die $self->error_message("Failed to run gatk IndelRealigner on normal");
    }

    return 1;
}

sub _run_varscan {
    my $self = shift;

    my $realigned_tumor_bam_file = $self->_realigned_tumor_bam_file;
    my $realigned_normal_bam_file = $self->_realigned_normal_bam_file;
    my $reference = $self->reference_fasta;
    my $output_indel = $self->varscan_indel_output;
    my $output_snp = $self->varscan_snp_output;
    my $varscan_params = $self->varscan_params;
    my $small_indel_list = $self->small_indel_output_bed;
    (my $small_indel_list_nobed = $small_indel_list) =~ s/\.padded1bp\.bed$/\.annotation_format/;
    my $final_output_file = $self->final_output_file;

    my $rv = Genome::Model::Tools::Varscan::Validation->execute(
        normal_bam => $realigned_normal_bam_file,
        tumor_bam => $realigned_tumor_bam_file,
        output_indel => $output_indel,
        output_snp => $output_snp,
        reference => $reference,
        varscan_params => $varscan_params,
    );
    die $self->error_message("Failed to run gmt varscan validation") unless $rv->result == 1;

    $rv = Genome::Model::Tools::Varscan::ProcessValidationIndels->execute(
        validation_indel_file => $output_indel,
        validation_snp_file => $output_snp,
        variants_file => $small_indel_list_nobed,
        output_file => $final_output_file,
    );
    die $self->error_message("Failed to run gmt varscan process-validation-indels") unless $rv->result == 1;

    return 1;
}

sub _resolve_output_directory {
    my $self = shift;
    my $output_dir;
    if ($self->build) {
        $output_dir = File::Spec->join($self->build->data_directory, 'indel_validation');
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
            || $self->large_indel_output_bed || $self->varscan_indel_output || $self->varscan_snp_output) {
            die $self->error_message("If a build is provided, you should not provide other params");
        }

        my $base_dir = $self->_resolve_output_directory();
        $self->final_output_file("$base_dir/final_output");
        $self->realigned_bam_file_directory("$base_dir/realigned_bams");
        $self->small_indel_output_bed("$base_dir/small_indels.bed");
        $self->large_indel_output_bed("$base_dir/large_indels.bed");
        $self->varscan_indel_output("$base_dir/varscan_indels");
        $self->varscan_snp_output("$base_dir/varscan_snps");
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
