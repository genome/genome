package Genome::InstrumentData::Command::List::Imported;

#rlong 2/22/10
# This command is picking up properties of InstrumentData
# as well as InstrumentData::Imported.
# This command is very slow. On average, it takes 4 minutes 
# to return a listing of mg.imported_instrument_data. I'm not sure why.

use strict;
use warnings;

use Genome;
use Command; 
use Data::Dumper;

class Genome::InstrumentData::Command::List::Imported {
    is => 'UR::Object::Command::List',
    has => [
        subject_class_name  => {
            is_constant => 1, 
            value => 'Genome::InstrumentData::Imported' 
        },
        show => { default_value => 'id,sample_name,sequencing_platform,import_format' },
    ],
    doc => 'list imported instrument data available for analysis',
};

1;

#$HeadURL: /gscpan/perl_modules/trunk/Genome/InstrumentData/Command/List.pm $
#$Id: /gscpan/perl_modules/trunk/Genome/InstrumentData/Command/List.pm 41086 2008-11-17T19:51:31.012449Z ebelter  $
