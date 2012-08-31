package Genome::Report::Dataset;
#:adukes is this general enough to live outside Report namespace?
#:ebelter @adukes - prolly

use strict;
use warnings;

use Genome;

use Carp 'confess';
use Data::Dumper 'Dumper';
use XML::LibXML;

class Genome::Report::Dataset {
    is => 'UR::Object',
    has => [
    name => { 
        is => 'Text',
        doc => 'Name of the dataset.',
    },
    row_name => {
        is => 'Text',
        doc => 'The name for each row of data.',
    },
    headers => {
        is => 'ARRAY',
        doc => 'Headers for the data.'
    },
    ],
    has_optional => [
    attributes => {
        is => 'HASH',
        default_value => {},
        doc => 'Attributes of the dataset.'
    },
    rows => {
        is => 'ARRAY',
        default_value => [],
        doc => 'Rows of data.'
    },
    ],
};

#< Rows >#
sub add_row {
    my ($self, $row) = @_;

    $self->_validate_aryref(
        name => 'row',
        value => $row,
        method => 'add row',
    )or return;
    
    return push @{$self->rows}, $row;
}

#< Attributes >#
sub get_attribute {
    my ($self, $name) = @_;

    confess $self->error_message("No name to get attribute") unless $name;
    
    return $self->attributes->{$name};
}

sub set_attribute {
    my ($self, $name, $value) = @_;

    $self->_validate_string_for_xml(
        name => 'name of attribute to set',
        value => $name,
        method => 'set attribute',
    ) or return;
    
    return $self->attributes->{$name} = $value;
}

#< Data Grab >#
sub get_row_values_for_header {
    my ($self, $header) = @_;

    unless ( @{$self->rows}) { 
        $self->error_message("No rows to get row vaalues from for header.");
        return;
    }

    unless ( $header ) {
        $self->error_message("No header given to get values from rows.");
        return;
    }

    my $pos;
    my $headers = $self->headers;
    for (my $i = 0; $i <= @{$self->headers}; $i++) {
        if ( $headers->[$i] eq $header ) {
            $pos = $i;
            last;
        }
    }
    unless ( defined $pos ) {
        $self->error_message("Header ($header) not found in headers.");
        return;
    }

    return map { $_->[$pos] } @{$self->rows};
}

#< Create >#
sub create {
    my $class = shift;

    my $self = $class->SUPER::create(@_)
        or return;

    # Names
    for my $attr (qw/ name row_name /) {
        $self->_validate_string_for_xml(
            name => join(' ', split('_', $attr)),
            value => $self->$attr,
            method => 'create',
        ) or return;
    }

    # Headers
    $self->_validate_aryref(
        name => 'headers',
        value => $self->headers,
        method => 'create',
    ) or return;
    for my $header ( @{$self->headers} ) {
        $self->_validate_string_for_xml(
            name => 'header',
            value => $header,
            method => 'create',
        ) or return;
    }

    # Rows
    if ( $self->rows ) {
        $self->_validate_aryref(
            name => 'rows',
            value => $self->rows,
            method => 'create',
        ) or return;
        for my $row ( @{$self->rows} ) {
            $self->_validate_aryref(
                name => 'row of rows',
                value => $row,
                method => 'create',
            ) or return;
        }
    }

    return $self;
}

sub create_from_xml_element {
    my ($class, $element) = @_;

    # Row Nodes
    my @row_nodes = grep { $_->nodeType == 1 } $element->findnodes('*');

    confess "No rows found in element to create dataset." unless @row_nodes;
    
    # Names
    my $name = $element->nodeName;
    my $row_name = $row_nodes[0]->nodeName;

    # Headers
    my $headers = [ map { 
        $_->nodeName 
    } grep {
        $_->nodeType == 1
    } $row_nodes[0]->findnodes('*') ];

    # Rows
    my @rows;
    for my $row ( @row_nodes ) {
        push @rows, [ map { $_->to_literal } grep { $_->nodeType == 1 } $row->getChildnodes ];
    }

    my %attributes;
    for my $attribute_node ( $element->attributes ) {
        $attributes{ $attribute_node->nodeName } = $attribute_node->getValue;
    }

    # Create
    return $class->create(
        name => $name,
        row_name => $row_name,
        headers => $headers,
        rows => \@rows,
        attributes => \%attributes,
    );
}

#< Validations >#
sub _validate_aryref { 
    my ($self, %params) = @_;

    # value => value of attr
    # name => name of attr
    # method => caller method

    unless ( $params{value} ) {
        $self->error_message(
            sprintf(
                '"%s" is/are required to "%s"',
                ucfirst($params{name}),
                $params{method},
            )
        );
        return;
    }

    my $ref = ref $params{value};
    unless ( $ref and $ref eq 'ARRAY' ) {
        $self->error_message(
            sprintf(
                '"%s" is/are required to be an array reference to "%s"',
                ucfirst($params{name}),
                $params{method},
            )
        );
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
        $self->error_message( 
            sprintf(
                '"%s" is/are required to "%s"',
                ucfirst($params{name}),
                $params{method}
            )
        );
        return;
    }

    if ( $params{value} =~ /\s/ ) {
        $self->error_message(
            sprintf(
                'Spaces were found in "%s" from method "%s", and are not allowed',
                $params{name},
                $params{method},
            )
        );
        return;
    }

    return 1;
}

#< XML - Create on the fly >#
sub to_xml_string {
    my $self = shift;
    
    my $node = $self->to_xml_element;
    unless ( $node ) {
        $self->error_message("Can't get xml element, which is necessarty to get xml string");
        return;
    }

    return $node->toString(1);
}

sub to_xml_element {
    my $self = shift;

    my $libxml = XML::LibXML->new();
    
    # Element
    my $element = XML::LibXML::Element->new( $self->name )
        or return;

    # Rows
    my $headers = $self->headers;
    for my $row ( @{$self->rows} ) {
        my $row_element = XML::LibXML::Element->new( $self->row_name )
            or return;
        $element->addChild($row_element)
            or return;
        for ( my $i = 0; $i < @$headers; $i++ ) {
            #print Dumper([$headers->[$i], $row->[$i]]);
            my $element = $row_element->addChild( XML::LibXML::Element->new($headers->[$i]) )
                or return;
            $element->appendTextNode($row->[$i])if(defined $row and defined $row->[$i]);
            #$row_element->addChild( $libxml->createAttribute($headers->[$i], $row->[$i]) )
            #    or return;
        }
    }

    # Add dataset attributes
    if ( my $attrs = $self->attributes ) {
        for my $attr ( keys %$attrs ) {
            $element->addChild( XML::LibXML::Attr->new($attr, $attrs->{$attr}) )
                or return;
        }
    }

    return $element;
}

#< SVS >#
sub to_separated_value_string {
    my ($self, %params) = @_;

    # separator
    my $separator = ( defined $params{separator} ) ? $params{separator} : ',';

    my $svs;
    
    # headers - default is to include them
    unless ( defined $params{include_headers} and $params{include_headers} == 0 ) {
        $svs = join($separator, @{$self->headers})."\n";
    }

    # rows
    for my $row ( @{$self->rows} ) {
        $svs .= join($separator, @$row)."\n";
    }

    return $svs;
}

1;

=pod

=head1 Name

Genome::Report::Dataset

=head1 Synopsis

This package is an interface to a set of data that has a name, attributes, headers and rows. Datasets give access to data creted directly or retreived from XML.  In a report generator, datasets can be created and added to a report via the 'add_dataset' method. It can also generate the xml string for the data, a separated value string of the data and added to an XML::LibXML doc. Reports can get their datasets, giving the same interface to the data.

=head1 Usage

 # Create (see properties below for explanations)
 my $dataset = Genome::Report::Dataset->create(
    name       => 'people',
    row_name   => 'person',
    headers    => [qw/ name age /],
    rows       => [ [qw/ Watson 80 /], [qw/ Crick NA /] ],
    attributes => { discovered => 'dna', },
 );
 
 # Rows
 $dataset->add_row([qw/ Franklin 75 /]); # Returns 1
 $dataset->get_row_values_for_header('age'); # Retruns (80, 'NA')
 
 # Attributes
 $dataset->get_attribute('discovered'); # Returns 'dna'
 $dataset->set_attribute('discovered', 'dna structure'); # Returns 'dna structure'
 
 # Convert
 $dataset->to_xml_element; # Returns an XML::LibXML::Element
 
 $dataset->to_xml_string; 
 # Returns the sting version of xml element above by calling toString(1) on it.
 <people discovered="dna">
  <person>
    <name>Watson</name>
    <age>80</age>
  </person>
  <person>
    <name>Crick</name>
    <age>NA</age>
  </person>
  <person>
    <name>Franklin</name>
    <age>70</age>
  </person>
 </people>

 $dataset->to_separated_value_string;
 # Returns:
 name,age
 Watson,80
 Crick,NA
 Franklin,75

=head1 Properties

B<Required>

=over 2

=item B<name>       The name of the dataset

=item B<row_name>   The name of the individual rows

=item B<headers>    The column headers - arrayref
 
=back

B<Optional>

=over

=item B<rows>       The rows of data - arrayref of arrayrefs

=item B<attributes> Meta data - hashref

=back

=head1 Public Methods

=head2 add_row
 
=over

=item I<Synopsis>   Add a row to the other rows of data.  No check is made that the number of elements in the row match the number of elements in the heaeder.

=item I<Arguments>  Arrayref.

=item I<Returns>    1 for success, 0 for failure.

=back

=head2 get_attribute

=over

=item I<Synopsis>   Get an attribute for the dataset.

=item I<Arguments>  Attribute's name.

=item I<Returns>    Attribute's value, if exists.

=back

=head2 set_attribute

=over

=item I<Synopsis>   Set an attribute on the dataset.

=item I<Arguments>  Attribute's name and value.

=item I<Returns>    Attribute's value.

=back

=head2 get_row_values_for_header

=over

=item I<Synopsis>   Get all the values for a particular header.

=item I<Arguments>  Header's name.

=item I<Returns>    Array of the row values.

=back

=head2 to_xml_element

=over

=item I<Synopsis>   Creates an XML element object the represents the dataset.  This element can then be added to an XML::LibXML element.

=item I<Arguments>  None.

=item I<Returns>    XML element object (XML::LibXml::Element) suitable for framing.

=back

=head2 to_xml_string

=over

=item I<Synopsis>   Creates and XML element as above, and then converts it to a string.

=item I<Arguments>  None.

=item I<Returns>    XML string.

=back

=head2 to_separated_value_string

=over

=item I<Synopsis>   Creates a string from the headers and rows.  Values are separated by the indicated separtor and rows are separated by a newline.  Does not include the dataset name or attributes.

=item I<Arguments>  Separator, defaults to comma (,). 

=item I<Returns>    Separated value string.

=back

=head1 Disclaimer

Copyright (C) 2009 Washington University Genome Sequencing Center

This module is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY or the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

=head1 Author(s)

B<Eddie Belter> I<ebelter@genome.wustl.edu>

=cut

#$HeadURL$
#$Id$
