package Genome::Model::Tools::DetectVariants2::Filter::PindelReadSupport;

use warnings;
use strict;

use Genome;

my %positions;

class Genome::Model::Tools::DetectVariants2::Filter::PindelReadSupport{
    is => 'Genome::Model::Tools::DetectVariants2::Filter',
    has => [
        min_variant_support => {
            is => 'String',
            is_optional => 1,
            default => '0',
            doc => 'Required number of variant-supporting reads. Note: Pindel doesn\'t actually report the indel if var-support < 3.',
        },
        var_to_ref_read_ratio => {
            is => 'String',
            is_optional => 1,
            default => '0.2',
            doc => 'This ratio determines what ratio of variant supporting reads to reference supporting reads to allow',
        },
        remove_single_stranded => {
            is => 'Boolean',
            is_optional => 1,
            default => 1,
            doc => 'enable this to filter out variants which have exclusively pos or neg strand supporting reads.',
        },
        sw_ratio => {
            is => 'String',
            is_optional => 1,
            is_input => 1,
            default => '0.25',
            doc => 'Throw out indels which have a normalized ratio of normal smith waterman reads to tumor smith waterman reads (nsw/(nsw+tsw)) at or below this amount.',
        },
        capture_data => {
            is => 'Boolean',
            doc => 'Set this to cause the read-support filter to dump reads for each indel to a temp file in order to avoid out of memory errors',
            is_optional => 1,
        },
    ],
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
    my $read_support_file = $self->_temp_staging_directory."/indels.hq.read_support.bed";
    my $hq_raw_file = $self->_temp_staging_directory."/indels.hq";
    my $big_output_file = $self->_temp_staging_directory."/indels.hq.big_output";
    my $lq_raw_file = $self->_temp_staging_directory."/indels.lq";
    my $output_file = $self->_temp_staging_directory."/indels.hq.bed";
    my $output_lq_file = $self->_temp_staging_directory."/indels.lq.bed";
    my $indel_file = $self->input_directory."/indels.hq.bed";
    my $capture_data = $self->capture_data;

    my $input_file = $self->input_directory."/indels.hq";
    my $ppr_cmd = Genome::Model::Tools::Pindel::ProcessPindelReads->create(
                    input_file => $input_file,
                    output_file => $output_file,
                    reference_build_id => $self->reference_build_id,
                    mode => 'read_support',
                    create_hq_raw_reads => 1,
                    aligned_reads_input => $self->aligned_reads_input,
                    control_aligned_reads_input => $self->control_aligned_reads_input,
                    lq_output_file => $output_lq_file,
                    hq_raw_output_file => $hq_raw_file,
                    lq_raw_output_file => $lq_raw_file,
                    read_support_output_file => $read_support_file,
                    capture_data => $capture_data,
                    big_output_file => $big_output_file, );
    unless($ppr_cmd->execute){
        die $self->error_message("Could not execute gmt pindel process-pindel-reads");
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
