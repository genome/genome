package Genome::InstrumentData::Command::Import::WorkFlow::ImportInstData;

use strict;
use warnings;

use Genome;

require File::Temp;

class Genome::InstrumentData::Command::Import::WorkFlow::ImportInstData {
    is => 'Command::V2',
    has_input => {
        work_flow_inputs => { is => 'Genome::InstrumentData::Command::Import::Inputs', },
    },
    has_output => {
        instrument_data => { is => 'Genome::InstrumentData', is_many => 1, },
    },
    has_optional_transient => {
        _working_directory => { is => 'Text', },
    },
};

sub shortcut {
    my $self = shift;
    my @instrument_data = $self->work_flow_inputs->instrument_data_for_original_data_path;
    return if not @instrument_data;
    $self->instrument_data(\@instrument_data);
    return 1;
}

sub execute {
    my $self = shift;
    $self->status_message('Import instrument data...');

    my $working_directory = $self->_resolve_working_directory;
    return if not $working_directory;

    my $space_available = $self->_verify_adequate_disk_space_is_available_for_source_files;
    return if not $space_available;

    my $dag = Genome::InstrumentData::Command::Import::WorkFlow::Builder->create(
        work_flow_inputs => $self->work_flow_inputs,
    )->build_workflow;
    return if not $dag;

    my $process = $self->work_flow_inputs->process;
    if ( $process ) {
        $dag->recursively_set_log_dir( $process->log_directory );
    }

    my $inputs = $self->work_flow_inputs->as_hashref;
    return if not $inputs;
    $inputs->{working_directory} = $self->_working_directory;

    my $success = $dag->execute_inline($inputs);
    die 'Run wf failed!' if not $success;

    $self->instrument_data($success->{instrument_data});

    $self->status_message('Import instrument data...done');
    return 1;
}

sub _resolve_working_directory {
    my $self = shift;

    my $tmp_dir = File::Temp::tempdir(CLEANUP => 1);
    if ( not $tmp_dir ) {
        $self->error_message('Failed to create tmp dir!');
        return;
    }

    return $self->_working_directory($tmp_dir);
}

sub _verify_adequate_disk_space_is_available_for_source_files {
    my $self = shift;
    my $helpers = Genome::InstrumentData::Command::Import::WorkFlow::Helpers->get;
    my $space_available = $self->work_flow_inputs->source_files->verify_adequate_disk_space_is_available_for_processing($self->_working_directory);
    return $space_available;
}

1;

