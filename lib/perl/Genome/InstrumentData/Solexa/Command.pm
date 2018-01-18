package Genome::InstrumentData::Solexa::Command;
use strict;
use warnings;

use Genome;

class Genome::InstrumentData::Solexa::Command {
    is => 'Command::Tree',
    english_name => 'genome instrument_data solexa command',
    has => [
        instrument_data => { is => 'Genome::InstrumentData::Solexa', id_by => 'instrument_data_id' },
        instrument_data_id => { is => 'Text', doc => 'identifies the instrument data by id' },
    ],
    doc => 'work with solexa instrument data',
};

use Genome::Command::Crud;
Genome::Command::Crud->init_sub_commands(
    target_class => 'Genome::InstrumentData::Solexa',
    target_name => 'instrument_data',
    target_name_pl => 'instrument data',
    create => {do_not_init => 1},
    delete => {do_not_init => 1},
    list => {do_not_init => 1},
    update => {only_if_null => 1, exclude => [qw/ library full_path ignored /]},
);

1;
