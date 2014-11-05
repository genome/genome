package Genome::VariantReporting::Command::Wrappers::SingleModel;

use strict;
use warnings;
use Genome;

class Genome::VariantReporting::Command::Wrappers::SingleModel {
    is => 'Genome::VariantReporting::Command::Wrappers::ModelPair',
    has => [
        plan_file_basename => {
            default_value => "cle_germline_report_TYPE.yaml",
        },
    ],
    has_calculated => [
        followup => {
            calculate => q| return undef|,
        },
    ],
};

sub get_aligned_bams {
    my $self = shift;
    my @aligned_bams;
    push @aligned_bams, $self->discovery->merged_alignment_result->id;
    return \@aligned_bams;
}

sub get_sample_and_bam_map {
    my $self = shift;

    return (
        $self->discovery->tumor_sample->name  => $self->discovery->tumor_bam,
    );
}

sub get_translations {
    my $self = shift;
    my %translations;
    $translations{normal} = $self->discovery->tumor_sample->name;
    if ($self->gold_sample_name) {
        $translations{gold} = $self->gold_sample_name;
    }
    return \%translations;
}

1;

