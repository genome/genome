# FIXME ebelter
#  remove
#
package Genome::Model::Command::List::VariantReviewLists;
use strict;
use warnings;

use Genome;

class Genome::Model::Command::List::VariantReviewLists {
    is => 'UR::Object::Command::List',
    has => [
        subject_class_name  => { is_constant => 1, value => 'Genome::VariantReviewList' },
        model               => { is_optional => 1 },
        show                => { default_value => 'id,author,rt_ticket,name' },
    ],
};

sub sub_command_sort_position { 10 }

1;
