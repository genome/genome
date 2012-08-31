package Genome::Model::Report::Table;

use strict;
use warnings;

use Genome;

use Data::Dumper 'Dumper';
use Lingua::EN::Inflect 'PL';

class Genome::Model::Report::Table {
    is => 'Genome::Report::Generator',
    has => [
        name => {
            is => 'Text',
            doc => 'Name of report.',
        },
        description => {
            is => 'Text',
            doc => 'Report description.',
        },
        headers => {
            is => 'ARRAY',
            doc => 'The headers (properties) of the data.',
        },
        rows => {
            is => 'ARRAY',
            doc => 'Rows of data to display.',
        },
    ],
    has_optional => [
        row_name => {
            is => 'Text',
            default_value => 'row',
        },
    ],
};

sub create { 
    my $class = shift;

    my $self = $class->SUPER::create(@_)
        or return;

    for my $property_name (qw/ name description /) {
        unless ( defined $self->$property_name ) {
            $self->error_message("Property ($property_name) is required.");
            return;
        }
    }

    for my $property_name (qw/ headers rows /) {
        my $value = $self->$property_name;
        unless ( $value ) {
            $self->error_message("Property ($property_name) is required.");
            return;
        }
        my $ref = ref($value);
        unless ( $ref and $ref eq 'ARRAY' ) {
            $self->error_message("Property ($property_name) is required to be an array reference.");
            return;
        }
    }

    return $self;
}

sub _add_to_report_xml {
    my $self = shift;

    # Add headers
    my $headers_node = $self->_xml->createElement('headers');
    $self->_main_node->addChild($headers_node);
    my @headers = @{$self->headers};
    for my $header ( @headers ) {
        my $node = $self->_xml->createElement('header');
        $headers_node->addChild($node);
        $header =~ s#[\-\_]#\ #g;
        $node->appendTextNode(  Genome::Utility::Text::capitalize_words($header) );
    }

    # Rows
    $self->_add_dataset(
        name => PL($self->row_name),
        row_name => $self->row_name,
        headers => $self->headers,
        rows => $self->rows,
    ) or return;

    return 1;
}

1;

#$HeadURL$
#$Id$
