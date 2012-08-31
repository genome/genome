package Genome::Nomenclature::View::Detail::Json;

use strict;
use warnings;
require UR;

use XML::Simple;
use JSON;

class Genome::Nomenclature::View::Detail::Json {
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

##        perspective => $self->perspective(),
#    my %view_args = (
#        subject_class_name => $self->subject_class_name,
#        perspective => 'detail',
#        toolkit => 'xml'
#    );
#    my $xml_view = $obj->create_view(%view_args);
#    my $xml = $xml_view->content();
#    my $hash = XMLin($xml);

    my $hash = {};
    my $nomenclature = $obj;

    $hash->{name} = $nomenclature->name;
    $hash->{id} = $nomenclature->id;
    $hash->{'empty_equivalent'} = $nomenclature->empty_equivalent;
    $hash->{fields} = [];

    my $ds = $UR::Context::current->resolve_data_sources_for_class_meta_and_rule(Genome::Nomenclature->__meta__);
    my $dbh = $ds->get_default_dbh;

    for my $field ($nomenclature->fields) {
        my $f = {}; 
        $f->{id} = $field->id;
        $f->{name} = $field->name;
        $f->{type} = $field->type;
        my @enums = $field->enumerated_values;
        $f->{enumerated_values} = [map {$_->value} @enums];
        $f->{enumerated_value_ids} = [map {$_->id} @enums];
        my @enum_value_use_counts;
        for my $e (@enums) {
           my ($use_count) = $dbh->selectrow_array("select count(*) from genome_subject_attribute where nomenclature=? and attribute_value=?",{},$field->id, $e->value);
           push @enum_value_use_counts, $use_count;
        }
        $f->{enumerated_value_use_counts} = \@enum_value_use_counts;

        my ($use_count) = $dbh->selectrow_array("select count(*) from genome_subject_attribute where nomenclature=?",{},$field->id);
        $f->{use_count} = $use_count;
    
        push @{$hash->{fields}}, $f;
    }

    return $self->_json->encode($hash);
}



1;
