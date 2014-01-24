package Genome::Model::Tools::Vcf::CreateCrossSampleVcf::FilterNonCalls;

use strict;
use warnings;

use Genome;

class Genome::Model::Tools::Vcf::CreateCrossSampleVcf::FilterNonCalls {
    is  => 'Genome::Command::Base',
    has_input => [
        build_clump => {
            is => 'Genome::Model::Tools::Vcf::CreateCrossSampleVcf::BuildClump',
        },
    ],
    has_optional_output => [
        filtered_vcf => {
            via => 'build_clump',
        },
    ],
};

sub execute {
    my $self = shift;

    my $cmd = Genome::Model::Tools::Vcf::FilterNonCalls->create(
        input_file => $self->build_clump->vcf_file,
        output_file => $self->filtered_vcf,
    );
    return $cmd->execute();
}

1;
