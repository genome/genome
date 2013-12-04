package Genome::Model::DeNovoAssembly::Command::ProcessInstrumentData;

use strict;
use warnings;

use Genome;

require File::Temp;

class Genome::Model::DeNovoAssembly::Command::ProcessInstrumentData {
    is => 'Command::V2',
    has_input => [
        build => { is => 'Genome::Model::Build::DeNovoAssembly',
            is_output => 1 },
        instrument_data => { is => 'Genome::InstrumentData' },
    ],
    has_constant => [
        lsf_resource => { is => 'Text',
            default_value => "-R 'select[type==LINUX64 && mem>32000 && gtmp>200] rusage[mem=32000:gtmp=200] span[hosts=1]' -M 32000000",
        },
    ],
};

sub shortcut {
    my $self = shift;

    my %params = $self->build->read_processor_params_for_instrument_data(
        $self->instrument_data);

    my $result = Genome::InstrumentData::SxResult->get_with_lock(%params);

    if($result) {
        $self->status_message('Using existing result ' .
            $result->__display_name__);
        return $self->link_result_to_build($result);
    }
    else {
        return;
    }
}

sub execute {
    my $self = shift;
    my $instrument_data = $self->instrument_data;

    $self->status_message(
        'Process instrument data ' . $instrument_data->__display_name__ .
        ' for '.$self->build->description);

    my %params = $self->build->read_processor_params_for_instrument_data(
        $instrument_data);

    my $result = Genome::InstrumentData::SxResult->get_or_create(%params);

    $self->link_result_to_build($result);

    $self->status_message('Process instrument data...OK');

    return 1;
}

sub link_result_to_build {
    my $self = shift;
    my $result = shift;

    $result->add_user(label => 'processed_reads', user => $self->build);

    foreach my $output_file ($result->read_processor_output_files) {
        Genome::Sys->create_symlink($result->output_dir . '/' . $output_file,
            $self->build->data_directory . '/' . $output_file);
    }

    Genome::Sys->create_symlink(
        $result->read_processor_output_metric_file,
        $self->build->data_directory . '/'
            . $result->read_processor_output_metric_file_base_name);

    Genome::Sys->create_symlink(
        $result->read_processor_input_metric_file,
        $self->build->data_directory.'/'
            . $result->read_processor_input_metric_file_base_name);

    return 1;
}

1;
