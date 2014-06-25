package Genome::Model::Build::View::DataTable::Json;

use strict;
use warnings;
require UR;

use XML::Simple;
use JSON;

class Genome::Model::Build::View::DataTable::Json {
    is => 'UR::Object::View::Default::Json',
    has_constant => [
        toolkit     => { value => 'json' },
    ],
    has_optional => [
        encode_options => { is => 'ARRAY', default_value => ['ascii', 'pretty', 'allow_nonref', 'canonical'], doc => 'Options to enable on the JSON object; see the documentation for the JSON Perl module' },
    ],
};


sub _generate_content {

    my ($self) = @_;
    my $build = $self->subject();

    if (!$build) {
        Carp::confess('This JSON view couldnt get the subject of the view. class='
                    , $self->subject_class_name
                    . ' id='
                    . $self->subject_id);
    }

    my $json = {};
    my @allocs = $build->associated_disk_allocations();

    for my $a (@allocs) {
        my $j = [];

        for my $property (qw/kilobytes_requested 
                             archivable 
                             is_archived 
                             disk_group_name 
                             absolute_path 
                             id/) {
            push @$j, $a->$property;
        }
        push @{$json->{'aaData'}}, $j;
    }

    $json->{'aoColumns'} = [
        {'mDataProp' => 'kilobytes'},
        {'mDataProp' => 'is_archived'},
        {'mDataProp' => 'archivable'},
        {'mDataProp' => 'disk_group'},
        {'mDataProp' => 'absolute_path'},
        {'mDataProp' => 'id'},
    ];
    return $self->_json->encode($json);
}

1;
