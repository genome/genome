package Genome::Model::Set::View::Metrics::Json;

use strict;
use warnings;
require UR;

use XML::Simple;
use JSON;

class Genome::Model::Set::View::Metrics::Json {
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

    my %all_metric_names;
    my %builds;
    my %model_metric_attrs;

    for my $model ($set->members()) {
        my $build = $model->last_complete_build;
        if ($build) {
            $builds{$model->id} = $build;
            for my $metric ($build->metrics)
            {
               $model_metric_attrs{$model->id}->{$metric->name} = $metric->value;
               $all_metric_names{$metric->name} = 1; 
            }
        }
    }
    
    my $ret = { };
    my @fields = map { s/\s+/_/g; $_; } sort keys %all_metric_names;

    $ret->{'aoColumns'} = [{"mDataProp" => 'name'},
                           {"mDataProp" => 'build id'},
                           map {{"mDataProp" => $_}} @fields];
    $ret->{'aaData'} = [];

    for my $model ($set->members()) {
        my $row = [$model->__display_name__(), (exists $builds{$model->id} ? sprintf("<a href='/view/genome/model/build/status.html?id=%s'>%s</a>", $builds{$model->id}->id, $builds{$model->id}->id) : 'N/A' )];
        for my $field (sort keys %all_metric_names) {
            my $val = (defined $model_metric_attrs{$model->id}->{$field} ? $model_metric_attrs{$model->id}->{$field} : 'N/A');
            push @$row, $val;
        }
        push @{$ret->{'aaData'}}, $row;
    }
    
    return $self->_json->encode($ret);
}



1;
