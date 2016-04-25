package Genome::Model::Build::SomaticInterface;

use strict;
use warnings;
use Genome;

class Genome::Model::Build::SomaticInterface {
    is => 'UR::Object',
    is_abstract => 1,
};

sub reference_sequence_build {
    my $self = shift;
    $self->fatal_message('Abstract: (reference_sequence_build) needs to be defined on class (%s)', $self->class);
}

sub individual {
    my $self = shift;
    $self->fatal_message('Abstract: (individual) needs to be defined on class (%s)', $self->class);
}

sub individual_common_name {
    my $self = shift;
    $self->fatal_message('Abstract: (individual_common_name) needs to be defined on class (%s)', $self->class);
}

sub snvs_annotated_variants_vcf_file {
    my $self = shift;
    $self->fatal_message('Abstract: (snvs_annotated_variants_vcf_file) needs to be defined on class (%s)', $self->class);
}

sub indels_detailed_variants_vcf_file {
    my $self = shift;
    return File::Spec->join($self->data_directory, 'variants', 'indels.detailed.vcf.gz');
}

sub indels_effects_file {
    my $self = shift;
    my $tier = shift;
    $self->fatal_message('Abstract: (indels_effects_file) needs to be defined on class (%s)', $self->class);
}

sub snvs_effects_file {
    my $self = shift;
    my $tier = shift;
    return File::Spec->join($self->data_directory, 'effects', "snvs.hq.novel.$tier.v2.bed");
}

1;
