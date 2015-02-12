package Genome::Model::Command::InstrumentData::AlignmentBams;

use strict;
use warnings;

use File::Spec;

use Genome;

class Genome::Model::Command::InstrumentData::AlignmentBams {
    is => 'Command::V2',
    doc => 'Lists the paths of instrument-data alignment BAMs for the provided build',
    has => [
        build => {
            is => 'Genome::Model::Build',
            shell_args_position => 1,
        },
    ],
    has_optional => [
        outdir => {
            is => 'FileSystemPath',
        },
    ],
};

sub help_detail {
    return "List the paths of the instrument data alignment BAMs for the provided build. If outdir is specified results are written to a file in outdir else results are written to STDOUT.";
}

sub execute {
    my $self = shift;
    my $build = $self->build;

    unless ($build->can('alignment_results_for_instrument_data')) {
        die $self->error_message('This command does not work for builds of type %s.', $build->type_name);
    }

    my @io_params = (
        headers => ['INSTRUMENT_DATA_ID', 'FLOW_CELL_ID', 'LANE', 'BAM_PATH', 'BAMQC_PATH'],
        separator => "\t",
        in_place_of_null_value => '-',
    );
    if($self->outdir) {
        push @io_params, output => File::Spec->join($self->outdir, $build->id.".instrumentdataalignmentbams.txt");
    }

    my $writer = Genome::Utility::IO::SeparatedValueWriter->create(@io_params);

    for my $instrument_data ($build->instrument_data) {
        my $data = {};

        $data->{INSTRUMENT_DATA_ID} = $instrument_data->id;
        $data->{FLOW_CELL_ID } = eval { $instrument_data->flow_cell_id };
        $data->{LANE} = eval { $instrument_data->lane };

        my ($alignment_result) = $build->alignment_results_for_instrument_data($instrument_data);

        $data->{BAM_PATH} = $alignment_result->output_dir . '/all_sequences.bam' if $alignment_result;
        $data->{BAMQC_PATH} = $self->_get_bamqc_path($build, $alignment_result);
        $writer->write_one($data);
    }

    return 1;
}

#fills a hash reference, lane_bamqcpath, key is the lane
sub get_lane_bamqc_path {
    my $self = shift;
    my $build = shift;

    my $lane_bamqc_path = {};
    for my $instrument_data ($build->instrument_data) {
        my ($alignment_result) = $build->alignment_results_for_instrument_data($instrument_data);
        #Get the latest bamqc result
        my $bamqc_path = $self->_get_bamqc_path($build, $alignment_result);
        $lane_bamqc_path->{$instrument_data->id} = $bamqc_path;
    }
    return $lane_bamqc_path;
}

sub _get_bamqc_path {
    my $self = shift;
    my $build = shift;
    my $alignment_result = shift;

    my @bamqc_results =  Genome::InstrumentData::AlignmentResult::Merged::BamQc->get(
        alignment_result_id => $alignment_result->id,
        builds => $build,
    );

    my $max = '0';
    my $bamqc_result;
    for(@bamqc_results) {
        my $earliest_time = $_->best_guess_date_numeric;
        if ($earliest_time > $max) {
            $max = $earliest_time;
            $bamqc_result = $_;
        }
    }
    my $bamqc_path = $bamqc_result ? $bamqc_result->output_dir : '-';
}

1;
