package Genome::Model::Tools::TechD::SummarizeCaptureBuilds;

use strict;
use warnings;

use Genome;

class Genome::Model::Tools::TechD::SummarizeCaptureBuilds {
    is => ['Command'],
    has => [
        build_ids => { is => 'Text', doc => 'A comma delimited list of build ids to summarize. Example: 12345,23456,34567,67890' },
        alignment_summary => { is => 'Text', doc => 'The output tsv file of consolidated alignment summaries.' },
        coverage_summary => { is => 'Text', doc => 'The output tsv file of consolidated coverage metrics.' },
        wingspan => { is_optional => 1, default_value => 0, doc => 'The wingspan value that was added to each region of interest'},
    ],
};
sub help_brief {
    'A tool to consolidate multiple build summaries into one file.'
}

sub help_synopsis {
    'Consolidate alignment and coverage summaries from multiple builds'
}

sub help_detail {
    return <<"EOS"
Capture builds generate specific metrics for on/off target alignment and duplication.  In addition, RefCov runs at multiple minimum depth filters to evaluate coverage of each region of interest.  General coverage metrics are then summarized at each minimum depth.  This tool will take multiple builds and consolidate all the metrics for each build subject into two output file.  One for alignment metrics and one for coverage metrics.
EOS
}

sub execute {
    my $self = shift;

    my @build_ids = split(',',$self->build_ids);
    unless (Genome::Model::Tools::TechD::ConvergeAlignmentSummaries->execute(
        build_ids => \@build_ids,
        output_file => $self->alignment_summary,
        wingspan => $self->wingspan,
    )->result) {
        die('Failed to generate alignment summary '. $self->alignment_summary ." for builds:\n\t". join("\n\t", @build_ids));
    }
    unless (Genome::Model::Tools::TechD::ConvergeCoverageStatsSummaries->execute(
        build_ids => \@build_ids,
        output_file => $self->coverage_summary,
        wingspan => $self->wingspan,
    )->result) {
        die('Failed to generate coverage summary '. $self->coverage_summary ." for builds:\n\t". join("\n\t", @build_ids));
    }
    return 1;
}


1;
