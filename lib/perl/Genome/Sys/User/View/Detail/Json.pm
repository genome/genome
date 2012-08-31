package Genome::Sys::User::View::Detail::Json;

use strict;
use warnings;
require UR;

use XML::Simple;
use JSON;

class Genome::Sys::User::View::Detail::Json {
    is => 'UR::Object::View::Default::Json',
    has_constant => [
        toolkit     => { value => 'json' },
    ],
    has_optional => [
        encode_options => { 
            is => 'ARRAY', 
            default_value => ['ascii', 'pretty', 'allow_nonref', 'canonical'], 
            doc => 'Options to enable on the JSON object; see the documentation for the JSON Perl module' },
    ],
};


sub _generate_content {

    my $self = shift;
    my $user = $self->subject();

    if (!$user) {
        Carp::confess('This JSON view couldnt get the subject of the view. class='
                    . $self->subject_class_name
                    . ' id='
                    . $self->subject_id);
    }

    my $h = {};
    my @projects = $user->projects();

    for my $p (@projects) {
#        my $c = $p->creator();
#        my $creator = $c ? $c->username : '-';
        push @{ $h->{'aaData'} }, [$p->parts_count, $p->name, $p->id];
    }

    $h->{'aoColumns'} = [
        {'mDataProp' => 'parts_count'},
        {'mDataProp' => 'name'},
        {'mDataProp' => 'id'},
    ];

    $h->{'sScrollX'} = '100%';

    return $self->_json->encode($h);
}



1;
