package Genome::Model::MutationalSignificance::Command::CalcCovg;

use strict;
use warnings;

use Genome;
use Data::Dumper;

class Genome::Model::MutationalSignificance::Command::CalcCovg {
    is => ['Command::V2'],
    has_input => [
        music_build => {
            is => 'Genome::Model::Build',
            doc => 'Build that is using this software result',
        },
        somatic_variation_build => {
            is => 'Genome::Model::Build::SomaticVariation',
            doc => 'Build with tumor-normal pair to be analyzed',
        },
        output_dir => {
            is => 'Text',
            is_output => 1,
            doc => "Directory where output files and subdirectories will be written",
        },
        reference_sequence => {
            is => 'Text',
            doc => "Path to reference sequence in FASTA format",
        },
        normal_min_depth => {
            is => 'Text',
            doc => "The minimum read depth to consider a Normal BAM base as covered",
        },
        tumor_min_depth => {
            is => 'Text',
            doc => "The minimum read depth to consider a Tumor BAM base as covered",
        },
        min_mapq => {
            is => 'Text',
            doc => "The minimum mapping quality of reads to consider towards read depth counts",
        },
        roi_file => {
            is => 'Text',
            doc => "Tab delimited list of ROIs [chr start stop gene_name] (See Description)",
        },
    ],
    has_output => [
        output_file => {
            is => 'Text',
            doc => 'Path to gene coverage file for this sample',
        },
    ],
};

sub help_synopsis {
    return <<HELP
This module calculates per-feature coverage for a sample, given the somatic variation build of that sample.
General usage:

 genome model mutational-significance calc-covg \\
    --somatic-variation-build build_id \\
    --output-dir output_dir/ \\
    --reference-sequence input_dir/all_sequences.fa \\
    --roi-file input_dir/all_coding_exons.tsv \\

HELP
}

sub help_detail {
    return <<HELP;
This module wraps the MuSiC calc-covg tool to calculate coverage for a single sample.
The base counts are taken from the tumor and normal bam files in the somatic variation build.
For more details on the calculation, see gmt music bmr calc-covg --help
HELP
}

sub shortcut {
    my $self = shift;

    my %params = $self->_collect_params;
    my $result = Genome::Model::Build::SomaticVariation::CalcCovgResult->get_with_lock(%params);

    if ($result) {
        $self->debug_message('Using existing result ' .
                    $result->__display_name__);
        return $self->_link_result_to_build($result);
    }
    else {
        return;
    }
}

sub execute {
    my $self = shift;

    $self->debug_message("CalcCovg for build ".$self->somatic_variation_build->id);

    my $result = Genome::Model::Build::SomaticVariation::CalcCovgResult->get_or_create($self->_collect_params);

    unless ($result) {
        $self->error_message("Failed to generate CalcCovgResult");
        return;
    }

    my $status = "CalcCovg done";
    $self->debug_message($status);
    return $self->_link_result_to_build($result);
}

sub _collect_params {
    my $self = shift;

    my %params = (
        roi_file => $self->roi_file,
        reference_sequence => $self->reference_sequence,
        somatic_variation_build => $self->somatic_variation_build,
        test_name => ($ENV{GENOME_SOFTWARE_RESULT_TEST_NAME} || undef),
    );

    $params{normal_min_depth} = $self->normal_min_depth;
    $params{tumor_min_depth} = $self->tumor_min_depth;
    $params{min_mapq} = $self->min_mapq;
    return %params;
}

sub _link_result_to_build {
    my $self = shift;
    my $result = shift;

    my $output_dir = $self->output_dir;
    my $sample_name = $self->somatic_variation_build->tumor_build->subject->extraction_label;

    unless (-d $output_dir) {
        Genome::Sys->create_directory($output_dir);
    }

    my $output_file = $output_dir."/".$sample_name.".covg";

    $self->output_file($output_file);

    Genome::Sys->create_symlink(join('/', $result->output_dir, $result->output_file), $output_file);
    $result->add_user(label => 'calc_covg', user => $self->music_build);
    return 1;
}

1;
