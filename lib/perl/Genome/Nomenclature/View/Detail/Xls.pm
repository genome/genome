package Genome::Nomenclature::View::Detail::Xls;

use strict;
use warnings;
require UR;

use XML::Simple;
use Spreadsheet::WriteExcel;

class Genome::Nomenclature::View::Detail::Xls {
    is => 'UR::Object::View::Default::Text',
    has_constant => [
        toolkit     => { value => 'xls' },
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

    open my $workbook_fh, '>', \my $workbook_content or die "Failed to open FH $!";
    my $writer = Spreadsheet::WriteExcel->new($workbook_fh);
    my $sheet = $writer->add_worksheet($obj->name);
    
    $sheet->write_string(0, 0, "Subject Name");
    for my $i (0...$#fields) {
        $sheet->write_string(0, $i+1, $fields[$i]);
    }

    $writer->close;
    return $workbook_content;
}



1;
