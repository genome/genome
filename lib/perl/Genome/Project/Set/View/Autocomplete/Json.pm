package Genome::Project::Set::View::Autocomplete::Json;


use strict;
use warnings;
require UR;

use XML::Simple;
use JSON;

class Genome::Project::Set::View::Autocomplete::Json {
    is => 'UR::Object::View::Default::Text',
    has_constant => [
        toolkit     => { value => 'json' },
    ],
    has_optional => [
        encode_options => { is => 'ARRAY', default_value => ['ascii', 'pretty', 'allow_nonref', 'canonical'], doc => 'Options to enable on the JSON object; see the documentation for the JSON Perl module' },
    ],
};

my $json;
sub _json {
    my ($self) = @_;
    return $json if defined $json;

    $json = JSON->new;
    foreach my $opt ( @{ $self->encode_options } ) {
        eval { $json = $json->$opt; };
        if ($@) {
            Carp::croak("Can't initialize JSON object for encoding.  Calling method $opt from encode_options died: $@");
        }
        if (!$json) {
            Carp::croak("Can't initialize JSON object for encoding.  Calling method $opt from encode_options returned false");
        }
    }
    return $json;
}

sub _generate_content {
    my $self = shift;

    my $set = $self->subject();

    if (!$set) {
        Carp::confess('This JSON view couldnt get the set.');
    }

    my @m = $set->members();
    my $hash = { names => [ map {$_->name} @m ]  };

    return $self->_json->encode($hash);
}


1;

