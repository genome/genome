package Genome::InstrumentData::Command::List::454;

#REVIEW fdu 11/20/2009
#If the list gets too long, filter and/or show property should be
#forced to use.

use strict;
use warnings;

use Genome;
use Command; 
use Data::Dumper;

class Genome::InstrumentData::Command::List::454 {
    is => 'UR::Object::Command::List',
    has => [
        subject_class_name  => {
            is_constant => 1, 
            value => 'Genome::InstrumentData::454' 
        },
        show => { default_value => 'id,run_name,region_number,index_sequence,sample_name' }
    ],
    doc => 'list 454 regions available for analysis',
};

sub _base_filter {
    'sequencing_platform=454'
}

1;

#$HeadURL: /gscpan/perl_modules/trunk/Genome/InstrumentData/Command/List.pm $
#$Id: /gscpan/perl_modules/trunk/Genome/InstrumentData/Command/List.pm 41086 2008-11-17T19:51:31.012449Z ebelter  $
