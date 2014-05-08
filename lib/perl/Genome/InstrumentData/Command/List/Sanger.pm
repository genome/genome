package Genome::InstrumentData::Command::List::Sanger;

use strict;
use warnings;

use Genome;

class Genome::InstrumentData::Command::List::Sanger {
    is => 'Genome::Object::Command::List',
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

