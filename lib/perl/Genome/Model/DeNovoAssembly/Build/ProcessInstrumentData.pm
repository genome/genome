package Genome::Model::DeNovoAssembly::Build::ProcessInstrumentData;

use strict;
use warnings;

use Genome;
use Genome::Model::DeNovoAssembly::SxReadProcessor;

class Genome::Model::DeNovoAssembly::Build::ProcessInstrumentData {
    is => 'Command::V2',
    has_input => [
        build => { 
            is => 'Genome::Model::Build::DeNovoAssembly',
            is_output => 1,
            doc => 'Build to process instrument data.',
        },
        instrument_data => { 
            is => 'Genome::InstrumentData',
            doc => 'Instrument data to be processed.',
        },
    ],
    has_output => [
        sx_result => {
            is => 'Genome::InstrumentData::SxResult',
            doc => 'Sx result from processing instrument data.',
        },
    ],
    has_constant => [
        lsf_resource => { 
            is => 'Text',
            default_value => "-R 'select[type==LINUX64 && mem>32000 && gtmp>200] rusage[mem=32000:gtmp=200] span[hosts=1]' -M 32000000",
        },
    ],
};

sub shortcut {
    my $self = shift;
    $self->debug_message('Process instrument data, attempting to shortcut...');

    my $sx_result_params = $self->_get_sx_result_params;
    return if not $sx_result_params;

    my $result = Genome::InstrumentData::SxResult->get_with_lock(%$sx_result_params);
    if ( $result ) {
        $self->sx_result($result);
        $self->debug_message('SX result found! '.$result->__display_name__);
        $self->debug_message('Shortcutting!');
        return 1;
    }
    else {
        $self->debug_message('No existing SX result found. Cannot shortcut!');
        return;
    }
}

sub execute {
    my $self = shift;
    $self->debug_message('Process instrument data...');

    my $sx_result_params = $self->_get_sx_result_params;
    return if not $sx_result_params;

    my $result = Genome::InstrumentData::SxResult->get_or_create(%$sx_result_params);
    if ( not $result ) {
        $self->error_message('Failed to create SX result!');
        return;
    }

    $self->sx_result($result);
    return 1;
}

sub _get_sx_result_params {
    my $self = shift;

    $self->debug_message('Build: '.$self->build->__display_name__);
    my $instrument_data = $self->instrument_data;
    $self->debug_message('Instrument data: '.$instrument_data->__display_name__);
    my $read_processor = $self->build->processing_profile->read_processor;
    $self->debug_message('Read processor: '.$read_processor);

    my $sx_processor = Genome::Model::DeNovoAssembly::SxReadProcessor->create(
        processor => $read_processor, # what if no processing is desired?
    );
    if ( not $sx_processor ) {
        $self->error_message('Failed to create SX read processor!');
        return;
    }
    $sx_processor->determine_processing( $self->build->instrument_data );

    my $sx_result_params = $sx_processor->sx_result_params_for_instrument_data($instrument_data);
    return if not $sx_result_params;
    $self->debug_message('SX reults params: '.Data::Dumper::Dumper($sx_result_params));

    return $sx_result_params;
}

1;

