package Genome::InstrumentData::Command::List::Sanger;

#REVIEW fdu 11/20/2009
#OK

use strict;
use warnings;

use Genome;
use Command; 
use Data::Dumper;

class Genome::InstrumentData::Command::List::Sanger {
    is => 'UR::Object::Command::List',
    has => [
        subject_class_name  => {
            is_constant => 1, 
            value => 'Genome::InstrumentData::Sanger' 
        },
        #show => { default_value => 'id,name,subject_name,processing_profile_name' },
    ],
    doc => 'list sanger/3730 runs (96-well) available for analysis',
};

sub create {

    my $class = shift;

    my $self = $class->SUPER::create(@_);

    unless (defined($self->filter())) {
        $self->error_message("You must provide a --filter option when listing sanger instrument data");  
        return;
    }
   
    return $self;

}

sub _base_filter {
    'sequencing_platform=sanger'
}

1;

#$HeadURL: /gscpan/perl_modules/trunk/Genome/InstrumentData/Command/List.pm $
#$Id: /gscpan/perl_modules/trunk/Genome/InstrumentData/Command/List.pm 41086 2008-11-17T19:51:31.012449Z ebelter  $
