package Genome::Model::Tools::DetectVariants2::Filter::PindelVafFilter;

use warnings;
use strict;

use Genome;

my %positions;

class Genome::Model::Tools::DetectVariants2::Filter::PindelVafFilter{
    is => 'Genome::Model::Tools::DetectVariants2::Filter',
    has => [
        variant_freq_cutoff => {
            is => 'Text',
            doc => "Variant Allele Frequency at or above which to pass on",
            default => "0.2",
            is_input => 1,
            is_optional => 1,
        },
        capture_data => {
            is => 'Boolean',
            doc => 'Set this to cause the read-support filter to dump reads for each indel to a temp file in order to avoid out of memory errors',
            default => 0,
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
    my $hq_raw_file = $self->_temp_staging_directory."/indels.hq";
    my $big_output_file = $self->_temp_staging_directory."/indels.hq.big_output";
    my $lq_raw_file = $self->_temp_staging_directory."/indels.lq";
    my $output_file = $self->_temp_staging_directory."/indels.hq.bed";
    my $output_lq_file = $self->_temp_staging_directory."/indels.lq.bed";
    my $indel_file = $self->input_directory."/indels.hq.bed";
    my $variant_freq_cutoff = $self->variant_freq_cutoff;
    my $capture_data = $self->capture_data;

    my $input_file = $self->input_directory."/indels.hq";
    my $ppr_cmd = Genome::Model::Tools::Pindel::ProcessPindelReads->create(
                    input_file => $input_file,
                    output_file => $output_file,
                    reference_build_id => $self->reference_build_id,
                    mode => 'vaf_filter',
                    create_hq_raw_reads => 1,
                    aligned_reads_input => $self->aligned_reads_input,
                    control_aligned_reads_input => $self->control_aligned_reads_input,
                    lq_output_file => $output_lq_file,
                    hq_raw_output_file => $hq_raw_file,
                    lq_raw_output_file => $lq_raw_file,
                    big_output_file => $big_output_file, 
                    variant_freq_cutoff => $variant_freq_cutoff,
                    capture_data => $capture_data, );
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
