package Genome::Model::RnaSeq::Command::SpliceJunctionSummary;

use strict;
use warnings;

use Genome;

my $DEFAULT_LSF_RESOURCE = "-R 'select[mem>=8000] rusage[mem=8000]' -M 8000000";

class Genome::Model::RnaSeq::Command::SpliceJunctionSummary {
    is => ['Command::V2'],
    has_input_output => [
        build_id => {},
    ],
    has => [
        build => { is => 'Genome::Model::Build', id_by => 'build_id', },
    ],
    has_param => [
        lsf_resource => { default_value => $DEFAULT_LSF_RESOURCE },
    ],
    has_optional_output => [
        splice_junction_summary_result_id => {
            is => 'Number',
            doc => 'ID of the result from running this step',
        },
    ],
    has_optional => [
        splice_junction_summary_result => {
            is => 'Genome::InstrumentData::AlignmentResult::Merged::SpliceJunctionSummary',
            id_by => 'splice_junction_summary_result_id',
            doc => 'The result from running this step',
        },
    ],
    doc => 'generate a summary of splice junctions based on pipeline results'
};

sub sub_command_category { 'pipeline steps' }

sub sub_command_sort_position { 6 }

sub help_synopsis {
    return <<EOS
    A splice-junction based utility to summarize alignments.
EOS
}

sub help_brief {
    return <<EOS
    A splice-junction based utility to summarize alignments.
EOS
}

sub help_detail {
    return <<EOS
This script will examine the alignment results directory and generate basic statistics on splice junction utilization.
EOS
}


sub shortcut {
    my $self = shift;

    my %params = $self->params_for_result;
    my $result = Genome::InstrumentData::AlignmentResult::Merged::SpliceJunctionSummary->get_with_lock(%params);

    if($result) {
        $self->status_message('Using existing result ' . $result->__display_name__);
        return $self->link_result_to_build($result);
    } else {
        return;
    }
}

sub execute {
    my $self = shift;
    my $build = $self->build;

    my %params = (
        $self->params_for_result,
        log_directory => $build->log_directory,
    );

    my $result = Genome::InstrumentData::AlignmentResult::Merged::SpliceJunctionSummary->get_or_create(%params);

    $self->link_result_to_build($result);


    return 1;
}


sub params_for_result {
    my $self = shift;
    my $build = $self->build;
    my $pp = $build->processing_profile;
    
    my $alignment_result = $build->alignment_result;
    my $bedtools_version = $pp->bedtools_version;

    unless ($alignment_result) {
        die $self->error_message('No alignment result found for build: '. $build->id);
    }

    my $result_users = Genome::SoftwareResult::User->user_hash_for_build($build);
    $result_users->{'splice-junction-summary'} = $build;

    return (
        alignment_result_id => $alignment_result->id,
        bedtools_version => $bedtools_version,
        test_name => ($ENV{GENOME_SOFTWARE_RESULT_TEST_NAME} || undef),
        users => $result_users,
    );
}

sub link_result_to_build {
    my $self = shift;
    my $result = shift;
    my $build = $self->build;

    my $link = $build->junctions_directory;
    Genome::Sys->create_symlink($result->output_dir, $link);
    
    $self->splice_junction_summary_result($result);

    return 1;
}







