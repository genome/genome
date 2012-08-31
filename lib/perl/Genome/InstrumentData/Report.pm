package Genome::InstrumentData::Report;

#REVIEW fdu 11/20/2009
#Remove Genome::Sys from base class list

use strict;
use warnings;

use Genome;

use Data::Dumper 'Dumper';

class Genome::InstrumentData::Report {
    is => ['Genome::Report::Generator','Genome::Sys'],
    has => [
        instrument_data => {
            is => 'Genome::InstrumentData',
            id_by => 'instrument_data_id'
        },
        instrument_data_id => {
            is => 'Integer', 
            doc=> 'Instrument Data id'
        },
    ],
};

sub create {
    my ($class, %params) = @_;

    unless ( $params{instrument_data_id} ) {
        $class->error_message("A id is required to create an instrument data report");
        return;
    }

    my $self = $class->SUPER::create(%params)
        or return;

    unless ( $self->instrument_data ) {
        $self->error_message( sprintf('Can\'t get a instrument data for instrument_data_id (%s)', $self->instrument_data_id) );
        return;
    }

    return $self;
}

#< Report Classes >#
sub get_generic_report_classes {
    my $type_name = shift;

    unless ( $type_name ) {
        Carp::confess("No instrument data sub type given\n"); 
        return;
    }

    return Genome::Sys::get_classes_in_subdirectory_that_isa(
        'Genome/InstrumentData/Report',
        'Genome::Report::Generator',
    );
}

sub get_report_classes_for_type_name {
    #print Dumper(\@_);
    my $type_name = shift;
    
    unless ( $type_name ) {
        Carp::confess("No instrument data sub type given\n"); 
    }

    return Genome::Sys::get_classes_in_subdirectory_that_isa(
        'Genome/InstrumentData/'.Genome::Utility::Text::string_to_camel_case($type_name).'/Report', 
        'Genome::Report::Generator',
    );
}

sub get_report_class_for_generic_report_name {
    my ($report_name) = @_;

    unless ( $report_name ) {
        Carp::confess("No report name given to get generic report class");
    }
    
    return 'Genome::InstrumentData::Report::'.Genome::Utility::Text::string_to_camel_case($report_name);
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
        'Genome::InstrumentData::%s::Report::%s', 
        Genome::Utility::Text::string_to_camel_case($type_name),
        Genome::Utility::Text::string_to_camel_case($report_name),
    );
}

#< Report Subclasses >#
sub get_generic_report_subclasses { 
    my $type_name = shift;

    my @classes = get_generic_report_classes($type_name)
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
    my $type_name = shift;

    my @subclasses = get_generic_report_subclasses($type_name)
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

    $self->_add_instrument_data_info
        or return;
    
    return $self->SUPER::generate_report;
}

sub _add_instrument_data_info {
    my $self = shift;

    my $instrument_data_node = $self->_xml->createElement('instrument-data-info')
        or return;
    $self->_main_node->addChild($instrument_data_node)
        or return;

    my %objects_attrs = (
        instrument_data => [
        qw/ id sequencing_platform run_name subset_name sample_name library_name  /,
        ],
    );
    for my $object ( keys %objects_attrs ) {
        for my $attr ( @{$objects_attrs{$object}} ) {
            my $value = $self->$object->$attr;
            $attr =~ s#\_#\-#g;
            my $element = $instrument_data_node->addChild( $self->_xml->createElement($attr) )
                or return;
            $element->appendTextNode( defined $value ? $value : '' );
        }
    }

    return 1;
}

1;

