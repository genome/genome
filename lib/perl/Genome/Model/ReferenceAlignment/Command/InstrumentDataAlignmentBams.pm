package Genome::Model::ReferenceAlignment::Command::InstrumentDataAlignmentBams;

use strict;
use warnings;

use Genome;

class Genome::Model::ReferenceAlignment::Command::InstrumentDataAlignmentBams {
    is => 'Genome::Command::Base',
    doc => "List the paths of the instrument data alignment BAMs for the provided build.",
    has => [
        build => {
            is => 'Genome::Model::Build::ReferenceAlignment',
            shell_args_position => 1,
        },
    ],
};


sub help_detail {
    return "List the paths of the instrument data alignment BAMs for the provided build.";
}


sub execute {
    my $self = shift;

    print join("\t", 'INSTRUMENT_DATA_ID', 'FLOW_CELL_ID', 'LANE', 'BAM_PATH') . "\n";
    for my $instrument_data ($self->build->instrument_data) {
        my $instrument_data_id = $instrument_data->id;
        my $flow_cell_id = eval { $instrument_data->flow_cell_id } || '-';
        my $lane = eval { $instrument_data->lane } || '-';
        my ($alignment_result) = $self->build->alignment_results_for_instrument_data($instrument_data);
        my $bam_path = ($alignment_result ? $alignment_result->output_dir . '/all_sequences.bam' : '-');
        print join("\t", $instrument_data_id, $flow_cell_id, $lane, $bam_path) . "\n";
    }

    return 1;
}

1;

