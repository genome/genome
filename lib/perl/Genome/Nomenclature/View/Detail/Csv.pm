package Genome::Nomenclature::View::Detail::Csv;

use strict;
use warnings;
require UR;

use XML::Simple;
use Spreadsheet::WriteExcel;

class Genome::Nomenclature::View::Detail::Csv {
    is => 'UR::Object::View::Default::Text',
    has_constant => [
        toolkit     => { value => 'csv' },
    ],
};


sub _generate_content {
    my $self = shift;


    my $obj = $self->subject();

    if (!$obj) {
        Carp::confess('This XLS view couldnt get the subject of the view. class='
                    , $self->subject_class_name
                    . ' id='
                    . $self->subject_id);
    }

    my @fields = sort map {$_->name} $obj->fields;

    unshift @fields, "Subject Name";
    return join ",", @fields;
}



1;
