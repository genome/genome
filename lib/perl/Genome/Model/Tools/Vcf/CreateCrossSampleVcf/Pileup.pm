package Genome::Model::Tools::Vcf::CreateCrossSampleVcf::Pileup;

use strict;
use warnings;

use Genome;
use Genome::File::Vcf::Reader;

class Genome::Model::Tools::Vcf::CreateCrossSampleVcf::Pileup {
    is  => 'Genome::Command::Base',
    has_input => [
        build_clump => {
            is => 'Genome::Model::Tools::Vcf::CreateCrossSampleVcf::BuildClump',
        },
        ref_fasta => {
            is => 'File',
        },
        use_bgzip => {
            is => 'Boolean',
        },
        samtools_version => {
            is => 'Text',
        },
        samtools_params => {
            is => 'Text',
        },
        region_file => {
            is => 'File',
        },
    ],
    has_optional_output => [
        pileup_file => {
            is => 'File',
        },
    ],
    has_param => [
        lsf_resource => {
            default => "-R 'select[model!=Opteron250 && type==LINUX64 && tmp>1000 && mem>16000] span[hosts=1] rusage[tmp=1000:mem=16000]' -M 1610612736",
        }
    ],
};

sub execute {
    my $self = shift;

    $self->pileup_file($self->build_clump->pileup_output_file);
    my $cmd = Genome::Model::Tools::Sam::Pileup->create(
        bam_file => $self->build_clump->bam_file,
        output_file => $self->build_clump->pileup_output_file,
        reference_sequence_path => $self->ref_fasta,
        use_bgzip => $self->use_bgzip,
        samtools_version => $self->samtools_version,
        samtools_params => $self->samtools_params,
        region_file => $self->region_file,
    );
    return $cmd->execute();
}

1;
