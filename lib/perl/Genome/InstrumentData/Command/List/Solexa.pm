package Genome::InstrumentData::Command::List::Solexa;

use strict;
use warnings;

use Genome;

class Genome::InstrumentData::Command::List::Solexa {
    is => 'Genome::Object::Command::List',
    has => [
        subject_class_name  => {
            is_constant => 1, 
            value => 'Genome::InstrumentData::Solexa' 
        },
        show => { default_value => 'id,flow_cell_id,lane,index_sequence,sample_name,library_name,read_length,is_paired_end,clusters,median_insert_size,sd_above_insert_size,target_region_set_name' },
    ],
    doc => 'list illumina/solexa lanes available for analysis',
};

sub _base_filter {
    'sequencing_platform=solexa'
}

1;

#$HeadURL: /gscpan/perl_modules/trunk/Genome/InstrumentData/Command/List.pm $
#$Id: /gscpan/perl_modules/trunk/Genome/InstrumentData/Command/List.pm 41086 2008-11-17T19:51:31.012449Z ebelter  $
