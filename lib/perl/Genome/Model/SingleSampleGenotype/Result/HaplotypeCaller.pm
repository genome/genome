package Genome::Model::SingleSampleGenotype::Result::HaplotypeCaller;

use strict;
use warnings;

use Genome;

use Genome::Utility::Text;

class Genome::Model::SingleSampleGenotype::Result::HaplotypeCaller {
    is => 'Genome::SoftwareResult::StageableSimple',
    has_input => {
        alignment_result => {
            is => 'Genome::InstrumentData::AlignedBamResult',
            doc => 'The alignments on which to run the haplotype caller',
        },
        intervals => {
            is_many => 1,
            is => 'Text',
            doc => 'regions to restrict this result',
        },
        gvcf_gq_bands => {
            is_many => 1,
            is => 'Number',
            doc => 'GQ thresholds for this result',
        },
        emit_reference_confidence => {
            is => 'Text',
            doc => 'ERC option for this result',
        },
        haplotype_caller_version => {
            is => 'Text',
            doc => 'Version of GATK to use',
        },
    },
    has => {
        vcf_file => {
            is => 'Text',
            is_calculated => 1,
            calculate_from => ['output_dir'],
            calculate => q{
                File::Spec->join($output_dir, $self->_vcf_filename);
            },
        },
    },
    doc => 'SoftwareResult wrapper for GATK HaplotypeCaller',
};


sub _vcf_filename {
    my $self = shift;

    my @intervals = $self->intervals;
    my $filename = (@intervals == 1)? Genome::Utility::Text::sanitize_string_for_filesystem($intervals[0]) : 'output';

    $filename .= '.g.vcf.gz';

    return $filename;
}

sub _run {
    my $self = shift;

    my $input_bam = $self->alignment_result->bam_path;
    my $reference = $self->alignment_result->reference_build->full_consensus_path('fa');
    my $output_dir = $self->temp_staging_directory;
    my $filename = $self->_vcf_filename;

    my $output_file = File::Spec->join($output_dir, $filename);


    my $cmd = Genome::Model::Tools::Gatk::HaplotypeCaller->create(
        input_bam => $input_bam,
        reference_fasta => $reference,
        output_vcf => $output_file,
        emit_reference_confidence => $self->emit_reference_confidence,
        intervals => [$self->intervals],
        gvcf_gq_bands => [sort { $a <=> $b } $self->gvcf_gq_bands],
        version => $self->haplotype_caller_version,
        java_interpreter => Genome::Sys->java_executable_path('1.7'),
    );
    $cmd->execute or die("Failed to run GATK HaplotypeCaller");

    return 1;
}

1;
