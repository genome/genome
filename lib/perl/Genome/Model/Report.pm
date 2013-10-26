package Genome::Model::Report;

use strict;
use warnings;

use Genome;

use Data::Dumper 'Dumper';

class Genome::Model::Report {
    #:adukes remove double inheritance, G:U:FileSystem is a pox, jpeck concurs
    is => ['Genome::Report::Generator','Genome::Sys'],
    has => [
        build => {
            is => 'Genome::Model::Build',
            id_by => 'build_id'
        },
        build_id => {
            is => 'Text',
            doc=> 'Build id'
        },
        model => {
            is => 'Genome::Model',
            via => 'build'
        },
        model_id => {
            via => 'model',
            to => 'id',
        },
        model_name => {
            via => 'model',
            to => 'name',
        },
    ],
};

sub create {
    my ($class, %params) = @_;

    unless ( $params{build_id} ) {
        $class->error_message("A build id is required to create a model report");
        return;
    }

    my $self = $class->SUPER::create(%params)
        or return;

    unless ( $self->build ) {
        $self->error_message( sprintf('Can\'t get a build for build_id (%s)', $self->build_id) );
        return;
    }

    return $self;
}

#< Report Classes >#
sub get_generic_report_classes {
    return Genome::Sys::get_classes_in_subdirectory_that_isa(
        'Genome/Model/Report',
        'Genome::Report::Generator',
    );
}

sub get_report_classes_for_type_name {
    my $type_name = shift;
    
    unless ( $type_name ) {
        Carp::confess("No model sub type given\n"); 
    }

    return Genome::Sys::get_classes_in_subdirectory_that_isa(
        'Genome/Model/'.Genome::Utility::Text::string_to_camel_case($type_name).'/Report', 
        'Genome::Report::Generator',
    );
}

sub get_report_class_for_generic_report_name {
    my ($report_name) = @_;

    unless ( $report_name ) {
        Carp::confess("No report name given to get generic report class");
    }
    
    return 'Genome::Model::Report::'.Genome::Utility::Text::string_to_camel_case($report_name);
}

sub get_report_class_for_type_name_and_report_name {
    my ($type_name, $report_name) = @_;

    unless ( $type_name ) {
        Carp::confess("No type name given to get report class");
    }

    unless ( $report_name ) {
        Carp::confess("No report name given to get report class");
    }

    return sprintf(
        'Genome::Model::%s::Report::%s', 
        Genome::Utility::Text::string_to_camel_case($type_name),
        Genome::Utility::Text::string_to_camel_case($report_name),
    );
}

#< Report Subclasses >#
sub get_generic_report_subclasses { 
    my @classes = get_generic_report_classes
        or return;

    return map { $_ =~ m#::([\w\d]+)$# } @classes;
}


sub get_report_subclasses_for_type_name {
    my $type_name = shift;

    my @classes = get_report_classes_for_type_name($type_name)
        or return;

    return map { $_ =~ m#::([\w\d]+)$# } @classes;
}

#< Report Names >#
sub get_generic_report_names { 
    my @subclasses = get_generic_report_subclasses
        or return;

    return map { Genome::Utility::Text::camel_case_to_string($_, ' ') } @subclasses;
}

sub get_report_names_for_type_name {
    my $type_name = shift;

    my @subclasses = get_report_subclasses_for_type_name($type_name)
        or return;

    return map { Genome::Utility::Text::camel_case_to_string($_, ' ') } @subclasses;
}

sub generate_report {
    my $self = shift;

    $self->_add_model_info
        or return;
    
    return $self->SUPER::generate_report;
}

sub _add_model_info {
    my $self = shift;

    my $build_node = $self->_xml->createElement('model-info')
        or return;
    $self->_main_node->addChild($build_node)
        or return;

    my %objects_attrs = (
        model => [
        qw/ id name type_name subject_name subject_type processing_profile_name /,
        $self->model->processing_profile->params_for_class 
        ],
        build => [qw/ build_id data_directory /],
    );
    for my $object ( keys %objects_attrs ) {
        for my $attr ( @{$objects_attrs{$object}} ) {
            if ($self->$object->can($attr)){
                my $value;
                my @values = $self->$object->$attr;
                my $property_meta = $self->$object->__meta__->property_meta_for_name($attr);
                if ($property_meta and $property_meta->is_many) {
                    $value = join(",",@values);
                }
                else {
                    $value = $values[0];
                }
                $attr =~ s#\_#\-#g;
                my $element = $build_node->addChild( $self->_xml->createElement($attr) )
                    or return;
                $element->appendTextNode( defined $value ? $value : '' );
            }
        }
    }

    return 1;
}

1;

#$HeadURL$
#$Id$
