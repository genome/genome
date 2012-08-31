package Genome::Task::View::Output::Json;

use strict;
use warnings;
require UR;

use XML::Simple;
use JSON;

class Genome::Task::View::Output::Json {
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
    $hash->{stdout} = $self->slurp_value($obj->stdout_pathname);
    $hash->{stderr} = $self->slurp_value($obj->stderr_pathname);
    $hash->{status} = $obj->status;
    $hash->{time_started} = $obj->time_started || "";
    $hash->{time_finished} = $obj->time_finished || "";

    return $self->_json->encode($hash);
}


sub slurp_value {
    my $self = shift;
    my $path = shift;
    my $val = "";
    
    if (-f $path) {
        $val = `cat $path`;
    }

    $val;
}


1;
