package Genome::Report::Generator;

use strict;
use warnings;

use Genome;

use Carp 'confess';
use Data::Dumper 'Dumper';
require Genome::Sys;
use IO::String;
use XML::LibXML;

class Genome::Report::Generator {
    is => 'UR::Object',
};

#< Create >#
sub create {
    my $class = shift;

    my $self = $class->SUPER::create(@_)
        or return;
    
    # XML doc
    $self->{_xml} = XML::LibXML->createDocument;
    unless ( $self->{_xml} ) {
        $self->error_message("Can't create XML object");
        $self->delete;
        return;
    }

    # main node
    $self->{_main_node} = $self->_xml->createElement('report')
        or Carp::confess('Create main node to XML');
    $self->_xml->addChild( $self->_main_node )
        or Carp::confess('Add main node to XML');
    $self->_xml->setDocumentElement( $self->_main_node );
    
    # meta node
    $self->{_meta_node} = $self->_xml->createElement('report-meta')
        or Carp::confess('Create meta node to XML');
    $self->_main_node->addChild( $self->_meta_node )
        or Carp::confess('Add meta node to main node');
    
    return $self;
}

#< Report Generation >#
sub generate_report {
    my $self = shift;

    # Subclass adds to report xml
    my $data;
    unless ( $data = $self->_add_to_report_xml ) {
        $self->status_message("No report data was generated.  This may be ok.  Errors, if any, would appear above.");
        return;
    }

    # Meta data
    $self->_add_report_meta_data
        or return;

    # Create report
    my $report = Genome::Report->create(
        xml => $self->_xml,
    );
    unless ( $report ) {
        $self->error_message("Can't create report.");
        return;
    }

    # DATA - BACKWARD COMPATIBILITY - THIS WILL BE REMOVED!
    if ( ref($data) ) {
        $self->warning_message("Generating a report w/ 'data' is deprecated.  Please store data as XML");
        $report->data($data);
    }

    return $report;
}

sub _add_report_meta_data {
    my $self = shift;

    # Basics
    for my $attr (qw/ name description date generator /) {
        my $value = $self->$attr;
        unless ( defined $value ) {
            $self->error_message("Report meta data attribute ($attr) is not defined.");
            return;
        }
        $self->_meta_node->addChild( $self->_xml->createElement($attr) )->appendTextNode($self->$attr);
    }

    # Generator params
    my %generation_params = $self->_get_params_for_generation;
    my $gen_params_node = $self->_xml->createElement('generator-params')
        or return;
    $self->_meta_node->addChild($gen_params_node)
        or return;
    for my $param ( keys %generation_params ) {
        for my $value ( @{$generation_params{$param}} ) {
            my $element = $gen_params_node->addChild( $self->_xml->createElement($param) )
                or return;
            $element->appendTextNode($value);
        }
    }

    return $self->_meta_node;
}

sub name {
    my ($subclass) = $_[0]->class =~ m#:?([\w\d]+)$#;
    confess "Can't get subclass from class: ".$_[0]->class unless $subclass;
    my $string =  Genome::Utility::Text::camel_case_to_string($subclass);
    return Genome::Utility::Text::capitalize_words($string);
}

sub description {
    return;
}

sub generator { 
    return $_[0]->class;
}

sub date { 
    return UR::Time->now; 
}

sub _get_params_for_generation {
    my $self = shift;

    my %params;
    for my $property ( $self->get_class_object->all_property_metas ) {
        my $property_name = $property->property_name;
        next if $property_name =~ /^_/;
        next if grep { $property_name eq $_ } (qw/ name description /);
        next if $property->via or $property->id_by;
        next unless $property->class_name->isa('Genome::Report::Generator');
        next if $property->class_name eq 'Genome::Report::Generator';
        #print Dumper($property_name);
        my $key = $property_name;
        $key =~ s#_#\-#g;
        next unless defined $self->$property_name;
        $params{$key} = [ $self->$property_name ];
    }

    return %params;
}

#< Helpers >#
sub _validate_aryref { 
    my ($self, %params) = @_;

    # value => value of attr
    # name => name of attr
    # method => caller method
    
    unless ( $params{value} ) {
        $self->error_message(ucfirst($params{name}).' are required for '.$params{method});
        return;
    }

    unless ( ref($params{value}) eq 'ARRAY' ) {
        $self->error_message(ucfirst($params{name}).' are required to be an array reference for '.$params{method});
        return;
    }

    return 1;
}

sub _validate_string_for_xml { 
    my ($self, %params) = @_;

    # value => value of attr
    # name => name of attr
    # method => caller method
    
    unless ( $params{value} ) {
        $self->error_message(ucfirst($params{name}).' is required for '.$params{method});
        return;
    }

    if ( $params{value} =~ /\s/ ) {
        $self->error_message(
            sprintf(
                'Spaces were found in %s from method %s, and are not allowed',
                $params{name},
                $params{method},
            )
        );
        return;
    }

    return 1;
}

#< XML >#
sub _xml {
    return $_[0]->{_xml};
}

sub _main_node {
    return $_[0]->{_main_node};
}

sub _meta_node {
    return $_[0]->{_meta_node};
}

sub _datasets_node {
    my $self = shift;

    unless ( $self->{_datasets_node} ) {
        $self->{_datasets_node} = $self->_xml->createElement('datasets')
            or return;
        $self->_main_node->addChild( $self->{_datasets_node} )
            or return;
    }

    return $self->{_datasets_node};
}

sub _add_dataset {
    my $self = shift;

    my $dataset;
    my $ref = ref($_[0]); # we gotta a dataset object
    if ( $ref ) {
        $dataset = $_[0];
    }
    else { # OLD way
        # separate the G:R:DS properties from the dataset xml attributes
        my %attrs = @_; 
        my %create_props = map { $_ => delete $attrs{$_} } (qw/ name row_name headers rows /);
        $dataset = Genome::Report::Dataset->create(
            %create_props,
            attributes => \%attrs,
        );
    }

    unless ( $dataset ) {
        $self->error_message();
        return;
    }
    
    my $element = $dataset->to_xml_element
        or return;
    $self->_datasets_node->addChild($element)
        or return;

    return $element;
}

sub _add_element_with_text_to_node {
    my ($self, $node, $name, $text) = @_;

    my $element = $node->addChild( $self->_xml->createElement($name) )
        or return;
    $element->appendTextNode($text)
        or return;

    return $element;
}

#< Report Template Files >#
sub get_xsl_file_for_html { 
    return $_[0]->_get_xsl_file_for_type('html');
}

sub _get_xsl_file_for_type { 
    my ($self, $type) = @_;

    my $module = $self->class;
    Carp::confess( # no xsl for base 
        "Can't get xsl file for base Genome::Report::Genertor class.  Get from subclass."
    ) if $module eq __PACKAGE__; 
    
    my $genome_dir = Genome->get_base_directory_name();
    my $inc_dir = substr($genome_dir, 0, -7);
    $module =~ s#::#/#g;

    return sprintf(
        '%s/%s.%s.xsl',
        $inc_dir,
        $module,
        $type
    );
}

1;

=pod

=head1 Name

Genome::Report::Generator

=head1 Synopsis

Base class for report generators.  Use this class as a base for yours...then implement a '_generate_data' method that returns a hashref.

=head1 Usage

 my $generator = Genome::Report::Generator->create(
    name => 'Happy', # required
    ..other params...
 );

 my $report = $generator->generate_report
    or die;

 $report->save('some_directory')
    or die;

=head1 Public Methods

=head2 generate_report

 my $report = $generator->generate_report
    or die;

=over

=item I<Synopsis>   Generates data and creates a Genome::Report

=item I<Arguments>  none

=item I<Returns>    Genome::Report

=back

=head1 Private Methods Implemented in Subclasses

=head2 _generate_data

=over

=item I<Synopsis>   Generates data and returns a hashref containing keys description, html (opt) and csv (opt)

=item I<Arguments>  none

=item I<Returns>    hashref 

=back

=head1 See Also

=head1 Disclaimer

Copyright (C) 2005 - 2008 Washington University Genome Sequencing Center

This module is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY or the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

=head1 Author(s)

B<Eddie Belter> I<ebelter@watson.wustl.edu>

=cut

#$HeadURL$
#$Id$
