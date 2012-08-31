package Genome::Model::Tools::DetectVariants2::Filter::PindelSomaticCalls;

use warnings;
use strict;

use Genome;

class Genome::Model::Tools::DetectVariants2::Filter::PindelSomaticCalls{
    is => 'Genome::Model::Tools::DetectVariants2::Filter',
    has_constant => [
        use_old_pindel => {
            is => 'Number',
            default => 1,
            doc => 'This will be updated when more than one version of pindel, 0.2,  is available',
        },
    ],
    has_param => [
        lsf_resource => {
            default => "-R 'span[hosts=1] rusage[mem=16000]' -M 1600000000",
        },
    ],
};

sub _variant_type { 'indels' };

sub _filter_variants {
    my $self = shift;

    my $output_file = $self->_temp_staging_directory."/indels.hq.bed";
    my $output_lq_file = $self->_temp_staging_directory."/indels.lq.bed";
    my $indel_file = $self->input_directory."/indels.hq.bed";

    $self->find_somatic_events($indel_file, $output_file, $output_lq_file);

    return 1;
}

sub find_somatic_events {
    my $self = shift;
    my $indel_file = shift;
    my $output_file = shift;
    my $output_lq_file = shift;
    my $raw_pindel_input = $self->detector_directory."/indels.hq";
    my $hq_raw_pindel_output = $self->_temp_staging_directory."/indels.hq";
    my $lq_raw_pindel_output = $self->_temp_staging_directory."/indels.lq";

    my $ppr_cmd = Genome::Model::Tools::Pindel::ProcessPindelReads->create(
                    input_file => $raw_pindel_input,
                    output_file => $output_file,
                    reference_build_id => $self->reference_build_id,
                    mode => 'somatic_filter',
                    sort_output => 1,
                    create_hq_raw_reads => 1,
                    hq_raw_output_file => $hq_raw_pindel_output,
                    lq_raw_output_file => $lq_raw_pindel_output,
                    lq_output_file => $output_lq_file, );
    unless($ppr_cmd->execute){
        die $self->error_message("Call to gmt pindel process-pindel-reads did not complete successfully");
    }
    unlink($output_file);
    return 1;
}

sub _check_file_counts {
    return 1;
}

sub _create_bed_file {
    my $self = shift;
    my $detector_file = shift;
    my $bed_file = shift;

    my $convert = Genome::Model::Tools::Bed::Convert::Indel::PindelToBed->create(
                    source => $detector_file,
                    output => $bed_file,
                    reference_build_id => $self->reference_build_id,
    );
    unless($convert->execute){
        die $self->error_message("Failed to convert detector to bed");
    }

    return 1;
}

1;
