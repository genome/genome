package Genome::Model::Tools::TechD::ConvergeCoverageStatsSummaries;

use strict;
use warnings;

use Genome;

#This is not a gmt, this should be a build command or a part of the build convergence process.
class Genome::Model::Tools::TechD::ConvergeCoverageStatsSummaries {
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

    my @cs;
    my @labels;
    for my $build_id ($self->build_ids) {
        my $build = Genome::Model::Build->get($build_id);
        unless ($build) {
            die('Failed to find build for id '. $build_id);
        }
        my $coverage_summary = $build->coverage_stats_summary_file($self->wingspan);
        unless ($coverage_summary && -s $coverage_summary) {
            die('Failed to find coverage summary for build '. $build->id);
        }
        push @cs, $coverage_summary;
        my $subject_name = $build->model->subject_name;
        push @labels, $subject_name;
    }
    unless (Genome::Model::Tools::BioSamtools::CompareCoverageStats->execute(
        input_files => \@cs,
        labels => \@labels,
        output_file => $self->output_file,
    )->result) {
        die('Failed to compare alignment summaries!');
    }
    return 1;
}

1;
