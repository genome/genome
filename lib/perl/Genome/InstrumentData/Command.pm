package Genome::InstrumentData::Command;
use strict;
use warnings;

use Genome;

class Genome::InstrumentData::Command {
    is => 'Command::Tree',
    english_name => 'genome instrument_data command',
    doc => 'work with instrument data',
};

use Genome::Command::Crud;
Genome::Command::Crud->init_sub_commands(
    target_class => 'Genome::InstrumentData',
    target_name => 'instrument_data',
    target_name_pl => 'instrument data',
    create => {do_not_init => 1},
    delete => {do_not_init => 1},
    list => {do_not_init => 1},
    update => {exclude => [qw/ 
        full_path import_format import_source_name 
        library run_name sequencing_platform subset_name
        /]},
);

1;
