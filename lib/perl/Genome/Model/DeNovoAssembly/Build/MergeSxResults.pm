package Genome::Model::DeNovoAssembly::Build::MergeSxResults;

use strict;
use warnings;

use Genome;
use Genome::Model::DeNovoAssembly::SxReadProcessor;

class Genome::Model::DeNovoAssembly::Build::MergeSxResults {
    is => 'Command::V2',
    has_input => [
        build => { 
            is => 'Genome::Model::Build::DeNovoAssembly',
            is_output => 1,
            doc => 'Build to process instrument data.',
        },
        sx_results => {
            is => 'Genome::InstrumentData::SxResult',
            is_many => 1,
            doc => 'Sx result from processing instrument data.',
        },
    ],
    has_output => [
        merged_sx_result => {
            is => 'Genome::InstrumentData::MergedSxResult',
            doc => 'Sx result from processing instrument data.',
        },
    ],
    has_constant => [
        lsf_resource => { 
            is => 'Text',
            default_value => "-R 'select[type==LINUX64 && tmp>200000] rusage[tmp=200000] span[hosts=1]'",
        },
    ],
};

sub shortcut {
    my $self = shift;
    $self->status_message('Merge SX results, attempting to shortcut...');

    my $sx_result_params = $self->_get_sx_result_params;
    return if not $sx_result_params;

    my $result = Genome::InstrumentData::MergedSxResult->get_with_lock(%$sx_result_params);
    if ( $result ) {
        $self->merged_sx_result($result);
        $self->status_message('Merged SX result found! '.$result->__display_name__);
        $self->status_message('Shortcutting!');
        return 1;
    }
    else {
        $self->status_message('No existing merged SX result found. Cannot shortcut!');
        return;
    }
}

sub execute {
    my $self = shift;
    $self->status_message('Merge SX results...');

    my $sx_result_params = $self->_get_sx_result_params;
    return if not $sx_result_params;

    my $result = Genome::InstrumentData::MergedSxResult->get_or_create(%$sx_result_params);
    if ( not $result ) {
        $self->error_message('Failed to create merged SX result!');
        return;
    }

    $self->merged_sx_result($result);
    $self->status_message('Merge SX results...done');
    return 1;
}

sub _get_sx_result_params {
    my $self = shift;

    $self->status_message('Build: '.$self->build->__display_name__);
    $self->status_message('Sx results: '.join(' ', map { $_->id } $self->sx_results));
    my @instrument_data = map { $_->instrument_data } $self->sx_results;
    $self->status_message('Instrument data: '.join(' ', map { $_->id } @instrument_data));
    my $read_processor = $self->build->processing_profile->read_processor;
    $self->status_message('Read processor: '.$read_processor);

    my $sx_processor = Genome::Model::DeNovoAssembly::SxReadProcessor->create(
        processor => $read_processor, # what if no processing is desired?
    );
    if ( not $sx_processor ) {
        $self->error_message('Failed to create SX read processor! '.$read_processor);
        return;
    }
    my $sx_result_params = $sx_processor->determine_sx_result_params_for_multiple_instrument_data(@instrument_data);
    return if not $sx_result_params;
    $self->status_message('SX result params: '.Data::Dumper::Dumper($sx_result_params));

    return $sx_result_params;
}

1;

