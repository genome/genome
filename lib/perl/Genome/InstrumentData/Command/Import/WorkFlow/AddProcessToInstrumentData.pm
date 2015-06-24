package Genome::InstrumentData::Command::Import::WorkFlow::AddProcessToInstrumentData;

use strict;
use warnings;

use Genome;

class Genome::InstrumentData::Command::Import::WorkFlow::AddProcessToInstrumentData { 
    is => 'Command::V2',
    has_input => {
        process => {
            is => 'Genome::InstrumentData::Command::Import::Process',
            doc => 'Process to add to newly created instrument data.',
        },
        instrument_data => {
            is => 'Genome::InstrumentData',
            is_many => 1,
            is_output => 1,
            doc => 'Newly created instrument data to add import process.',
        }
    },
};

sub execute {
    my $self = shift;
    $self->debug_message('Add process to instrument data...');

    $self->debug_message('Process: '. $self->process->id);
    for my $instrument_data ( $self->instrument_data ) {
        $self->debug_message('Instrument Data: '. $instrument_data->id);
        $instrument_data->add_attribute(
            attribute_label => 'process_id',
#my $import_file = File::Spec->join($data_dir, 'info.tsv');
            attribute_value => $self->process->id,
        );
    }

    $self->debug_message('Add process to instrument data...OK');
    return 1;
}

1;

