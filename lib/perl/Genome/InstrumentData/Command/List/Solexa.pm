package Genome::InstrumentData::Command::List::Solexa;

#REVIEW fdu 11/20/2009
#Filter property should be forced to use like that in
#G::I::C::L::Sanger to avoid long endless useless stdout list,
#especially with growing amount of solexa data stored. User should
#never be allowed to run this command without --filter option

use strict;
use warnings;

use Genome;
use Command; 
use Data::Dumper;

class Genome::InstrumentData::Command::List::Solexa {
    is => 'UR::Object::Command::List',
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
