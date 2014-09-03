package Genome::VariantReporting::Framework::Command::Wrappers::SingleModel;

use strict;
use warnings;
use Genome;

class Genome::VariantReporting::Framework::Command::Wrappers::SingleModel {
    is => 'Genome::VariantReporting::Framework::Command::Wrappers::ModelPair',
    has => [
        plan_file_basename => {
            default_value => "cle_germline_report_TYPE.yaml",
        },
    ],
    has_calculated => [
        validation => {
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

sub get_translations {
    my $self = shift;
    my %translations;
    $translations{d30_normal} = $self->discovery->tumor_sample->name;
    return \%translations;
}

sub report_names {
    return qw(cle_germline_simple_report cle_germline_full_report);
}
1;

