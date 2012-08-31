# FIXME ebelter
#  remove
#
package Genome::Model::Command::List::Events;

use strict;
use warnings;

use Genome;

use Command; 
use Data::Dumper;

class Genome::Model::Command::List::Events {
    is => 'UR::Object::Command::List',
    has => [
    subject_class_name  => {
         is_constant => 1, 
        value => 'Genome::Model::Event' 
    }, 
    ],
};

sub sub_command_sort_position { 12 }

1;

#$HeadURL$
#$Id$
