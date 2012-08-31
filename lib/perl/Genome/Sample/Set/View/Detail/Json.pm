package Genome::Sample::Set::View::Detail::Json;

use strict;
use warnings;
require UR;

use XML::Simple;
use JSON;

class Genome::Sample::Set::View::Detail::Json {
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

    my $set = $self->subject();

    if (!$set) {
        Carp::confess('This JSON view couldnt get the subject of the view. class='
                    . $self->subject_class_name
                    . ' id='
                    . $self->subject_id);
    }

    my $h = {};

    my %all_nomenclature_fields;
    my %sample_attrs;

    for my $sample ($set->members()) {
        for my $attr ($sample->attributes()) {
           $sample_attrs{$sample->id}->{$attr->nomenclature} = $attr->attribute_value;
           $all_nomenclature_fields{$attr->nomenclature_id}->{$attr->nomenclature} = 1; 
        }
    }
    
    my $ret = { };
    my @fields;
    for my $nomenclature_id (keys %all_nomenclature_fields) {
        my @ind_field_ids = keys %{$all_nomenclature_fields{$nomenclature_id}};
        my @ind_fields = sort {$a->name <=> $b->name} Genome::Nomenclature::Field->get(id=>\@ind_field_ids);
        push @fields, @ind_fields;
    }

    #my @fields = Genome::Nomenclature::Field->get(id=>[keys %all_nomenclature_fields]);
    $ret->{'aoColumns'} = [{"mDataProp" => 'name'},
                           {"mDataProp" => 'id'},
                           map {{"mDataProp" => $_->name}} @fields];
    $ret->{'aaData'} = [];

    for my $sample ($set->members()) {
        my $row = [$sample->name(), $sample->id];
        for my $field (@fields) {
            my $attr = Genome::SubjectAttribute->get(nomenclature=>$field->id, subject=>$sample);
            my $val = (defined $attr ? $attr->attribute_value : '-');
            push @$row, $val;
        }
        push @{$ret->{'aaData'}}, $row;
    }
    
=cut
    $h->{'aoColumns'} = [
        { "mDataProp" => "engine" },
        { "mDataProp" => "browser" },
        { "mDataProp" => "platform" },
        { "mDataProp" => "details.0" },
        { "mDataProp" => "details.1" }
    ];

    for my $sample ($set->members()) {
        my @attr;

        for my $a ($sample->attributes()) {

            push @attr, 
            push @{$h->{'aaData'}}, [
                $sample->name(),
                $sample->id(),
                $a->nomenclature_name(),
                $a->nomenclature_field_name(),
                $a->attribute_value()
            ];
        }
    }

=cut
    return $self->_json->encode($ret);
}



1;
