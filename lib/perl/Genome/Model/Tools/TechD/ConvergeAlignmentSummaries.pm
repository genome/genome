package Genome::Model::Tools::TechD::ConvergeAlignmentSummaries;

use strict;
use warnings;

use Genome;

#This is not a gmt, this should be a build command or a part of the build convergence process.
class Genome::Model::Tools::TechD::ConvergeAlignmentSummaries {
    is => ['Command'],
    has => [
        build_ids => {
            is_many => 1,
        },
        output_file => {
        },
        wingspan => {
            is_optional => 1,
            default_value => 0,
        }
    ],
};

sub execute {
    my $self = shift;

    my @as;
    my @labels;
    for my $build_id ($self->build_ids) {
        my $build = Genome::Model::Build->get($build_id);
        unless ($build) {
            die('Failed to find build for id '. $build_id);
        }
        # TODO: Eventually use the alignment_summary_hash_ref since percentages are calculated by that method
        my ($alignment_summary) = $build->alignment_summary_file($self->wingspan);
        unless ($alignment_summary && -s $alignment_summary) {
            die('Failed to find alignment summary for build '. $build->id);
        }
        push @as, $alignment_summary;
        my $subject_name = $build->model->subject_name;
        push @labels, $subject_name;
    }
    unless (Genome::Model::Tools::BioSamtools::CompareAlignmentSummaries->execute(
        input_files => \@as,
        labels => \@labels,
        output_file => $self->output_file,
    )->result) {
        die('Failed to compare alignment summaries!');
    }
    return 1;
}

1;
