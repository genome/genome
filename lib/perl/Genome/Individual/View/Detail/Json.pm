package Genome::Individual::View::Detail::Json;

use strict;
use warnings;
require UR;

use XML::Simple;
use JSON;

class Genome::Individual::View::Detail::Json {
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

    my $obj = $self->subject();

    if (!$obj) {
        Carp::confess('This JSON view couldnt get the subject of the view. class='
                    , $self->subject_class_name
                    . ' id='
                    . $self->subject_id);
    }

    my $hash = {};
    my @attributes;
    for my $a ($obj->attributes()) {
        push @attributes, [
            $a->nomenclature_name() || '[default]',
            $a->nomenclature_id() || 'null',
            $a->attribute_label() || '-',
            $a->attribute_value() || '-'
        ];
    }
    $hash->{'aaData'} = \@attributes;
    $hash->{'aoColumns'} = ['foo','bar'];

    return $self->_json->encode($hash);
}



1;
