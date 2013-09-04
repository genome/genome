package Genome::Model::ReferenceAlignment::Command::InstrumentDataAlignmentBams;

use strict;
use warnings;

use Genome;
use Date::Manip;

class Genome::Model::ReferenceAlignment::Command::InstrumentDataAlignmentBams {
    is => 'Genome::Command::Base',
    doc => "List the paths of the instrument data alignment BAMs for the provided build.",
    has => [
        build_id => {
            is => 'Number',
            shell_args_position => 1,
        },
    ],
};


sub help_detail {
    return "List the paths of the instrument data alignment BAMs for the provided build.";
}


sub execute {
    my $self  = shift;
    my $build = Genome::Model::Build->get($self->build_id);

    die $self->error_message('Please provide valid build id') unless $build;
    die $self->error_message('The provided build '.$build->id. ' is not reference alignment build.')
        unless $build->model->type_name eq 'reference alignment';

    print join("\t", 'INSTRUMENT_DATA_ID', 'FLOW_CELL_ID', 'LANE', 'BAM_PATH', 'BAMQC_PATH') . "\n";
    for my $instrument_data ($build->instrument_data) {
        my $instrument_data_id = $instrument_data->id;
        my $flow_cell_id = eval { $instrument_data->flow_cell_id } || '-';
        my $lane = eval { $instrument_data->lane } || '-';
        my ($alignment_result) = $build->alignment_results_for_instrument_data($instrument_data);
        my $bam_path = $alignment_result ? $alignment_result->output_dir . '/all_sequences.bam' : '-';
        #Get the latest bamqc result
        my @bamqc_results =  Genome::InstrumentData::AlignmentResult::Merged::BamQc->get(
            alignment_result_id => $alignment_result->id
        );

        my $max = '0';
        my $bamqc_result;
        for(@bamqc_results) {
            my $earliest_time = UnixDate($_->best_guess_date, "%s");
            if ($earliest_time > $max) {
                $max = $earliest_time;
                $bamqc_result = $_;
            }
        }

        my $bamqc_path = $bamqc_result ? $bamqc_result->output_dir : '-';
        print join("\t", $instrument_data_id, $flow_cell_id, $lane, $bam_path, $bamqc_path) . "\n";
    }

    return 1;
}

1;

