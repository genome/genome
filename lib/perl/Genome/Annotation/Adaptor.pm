package Genome::Annotation::Adaptor;

use strict;
use warnings;
use Genome;

class Genome::Annotation::Adaptor {
    is => 'Command::V2',
    has_input => [
        build => {
            is => 'Genome::Model::Build::RunsDV2',
        },
        variant_type => {
            is => 'Text',
            valid_values => ['snvs', 'indels'],
        }
    ],
    has_output => [
        bam_results => {
            is_many => 1,
            is => 'Genome::InstrumentData::AlignedBamResult',
        },
        vcf_result => {
            is => 'Genome::Model::Tools::DetectVariants2::Result::Vcf',
        },
    ],
};

sub execute {
    my $self = shift;
    $self->bam_results($self->resolve_bam_results);
    $self->vcf_result($self->resolve_vcf_result);
    return 1;
}

sub resolve_bam_results {
    my $self = shift;

    if ($self->build->isa('Genome::Model::Build::SomaticVariation')) {
        return $self->_resolve_bam_results_variation;
    } elsif ($self->build->isa('Genome::Model::Build::SomaticValidation')) {
        return $self->_resolve_bam_results_validation;
    } else {
        die "This adaptor can only work on SomaticValidation or SomaticVariation type builds";
    }
}

sub _resolve_bam_results_variation {
    my $self = shift;
    my @bam_results;
    for my $type qw(normal_build tumor_build) {
        push @bam_results, $self->build->$type->merged_alignment_result;
    }
    return \@bam_results;
}

sub _resolve_bam_results_validation {
    my $self = shift;
    return [ $self->build->control_merged_alignment_result, $self->build->merged_alignment_result ];
}

sub resolve_vcf_result {
    my $self = shift;

    my $accessor = sprintf('get_detailed_%s_vcf_result', $self->variant_type);
    my $result = eval {$self->build->$accessor};
    my $error = $@;
    if ($error) {
        $self->debug_message("No %s result found on build %s:\n%s",
            $self->variant_type, $self->build->id, $error);
    }
    return $result;
}

1;
