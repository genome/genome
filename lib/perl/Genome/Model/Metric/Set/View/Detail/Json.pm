package Genome::Model::Metric::Set::View::Detail::Json;

use strict;
use warnings;
require UR;

use XML::Simple;
use JSON;

class Genome::Model::Metric::Set::View::Detail::Json {
    is => 'UR::Object::View::Default::Json',
    has_constant => [
        toolkit     => { value => 'json' },
    ],
    has_optional => [
        encode_options => { is => 'ARRAY', default_value => ['ascii', 'pretty', 'allow_nonref', 'canonical'], doc => 'Options to enable on the JSON object; see the documentation for the JSON Perl module' },
    ],
};

# using this stuff to fill a datatable via ajax call

sub _generate_content {

    my ($self) = @_;
    my $set = $self->subject();

    if (!$set) {
        Carp::confess('This JSON view couldnt get the subject of the view. class='
                    , $self->subject_class_name
                    . ' id='
                    . $self->subject_id);
    }

    my @rows;
    for my $m ($set->members) {
        push @rows, [$m->build_id, $m->name, $m->value];
    }

    return $self->_json->encode({ 'aaData' => \@rows });
}



1;
