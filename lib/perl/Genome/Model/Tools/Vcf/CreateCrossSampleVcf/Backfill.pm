package Genome::Model::Tools::Vcf::CreateCrossSampleVcf::Backfill;

use strict;
use warnings;

use Genome;

class Genome::Model::Tools::Vcf::CreateCrossSampleVcf::Backfill {
    is  => 'Genome::Command::Base',
    has_input => [
        build_clump => {
            is => 'Genome::Model::Tools::Vcf::CreateCrossSampleVcf::BuildClump',
        },
        pileup_file => {
            is => 'Text',
        },
        use_bgzip => {
            is => 'Boolean',
        },
        region_file => {
            is => 'File',
        },
    ],
    has_optional_output => {
        backfilled_vcf => {
            is => 'File',
        },
    },
};

sub execute {
    my $self = shift;

    $self->backfilled_vcf($self->build_clump->backfilled_vcf);
    my $cmd = Genome::Model::Tools::Vcf::Backfill->create(
        output_file => $self->build_clump->backfilled_vcf,
        merged_positions_file => $self->region_file,
        pileup_file => $self->pileup_file,
        vcf_file => $self->build_clump->filtered_vcf,
        use_bgzip => $self->use_bgzip,
    );
    return $cmd->execute();
}

1;
