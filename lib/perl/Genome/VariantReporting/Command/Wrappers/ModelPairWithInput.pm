package Genome::VariantReporting::Command::Wrappers::ModelPairWithInput;

use strict;
use warnings;
use Genome;

class Genome::VariantReporting::Command::Wrappers::ModelPairWithInput {
    is => 'Genome::VariantReporting::Command::Wrappers::ModelPair',
    has => [
        other_snvs_vcf_input => {is => 'Path'},
        other_indels_vcf_input => {is => 'Path'},
    ],
};

sub input_vcf {
    my ($self, $variant_type) = @_;
    my $accessor = join("_", "other", $variant_type, "vcf_input");
    return $self->$accessor;
}

1;

