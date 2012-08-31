package Genome::Task::Set::View::DataTable::Json;

use strict;
use warnings;

use Genome;

class Genome::Task::Set::View::DataTable::Json {
    is => 'UR::Object::View::Default::Json',
    has_constant => [
        toolkit     => { value => 'json' },
    ],
    has_optional => [
        encode_options => { is => 'ARRAY', default_value => ['ascii', 'pretty', 'allow_nonref', 'canonical'], doc => 'Options to enable on the JSON object; see the documentation for the JSON Perl module' },
    ],
};

sub _generate_content {
    my $self = shift;

    my @set; 

    for my $obj ($self->subject->members) {
        push @set, [$obj->__display_name__, $obj->status, $obj->user_id, $obj->time_submitted, $obj->time_started, $obj->time_finished ]
    }
    
    return $self->_json->encode({aaData => \@set, 
            aoColumns => [
                    {mDataProp => 'command class',},
                    {mDataProp => 'status',},
                    {mDataProp => 'user id',},
                    {mDataProp => 'time submitted',},
                    {mDataProp => 'time started',},
                    {mDataProp => 'time finished',},
            ],
            wutgiTaskIds => [map {$_->id} $self->subject->members],
     });

}


1;
