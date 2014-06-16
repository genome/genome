package Genome::VariantReporting::GetInitialVcfResult;

use strict;
use warnings FATAL => 'all';
use Genome;

class Genome::VariantReporting::GetInitialVcfResult {
    is => 'Command::V2',
    has_input => [
        build_id => {
            is => 'Text',
        },
        variant_type => {
            is => 'Text',
            valid_values => ['snvs', 'indels'],
        },
    ],
    has_output => [
        output_result => {
            is => 'Genome::SoftwareResult',
        },
    ],
};

sub shortcut {
    my $self = shift;
    return $self->execute();
}

sub execute {
    my $self = shift;

    $self->output_result($self->build->get_detailed_vcf_result($self->variant_type));
    return 1;
}

sub build {
    my $self = shift;

    my $build = Genome::Model::Build->get($self->build_id);
    if ($build) {
        return $build;
    } else {
        die $self->error_message("Couldn't find a build for id (%s)",
            $self->build_id);
    }
}


1;
