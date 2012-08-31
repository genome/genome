package Genome::View::Static::Html;

use strict;
use warnings;

use Genome;

class Genome::View::Static::Html {
    is => 'UR::Object::View::Default::Xsl',
    is_abstract => 1,
    has => [
        perspective => {
            is_constant => 1,
            value => 'static',
        },
        desired_perspective => {
            default_value => 'default',
        },
        transform => {
            is_constant => 1,
            value => 1,
        },
        _xml_doc => {
            is => 'XML::LibXML::Document',
            doc => 'The LibXML document used to create the content for this view',
            is_transient => 1
        },
    ],
    doc => 'Classes that want to produce HTML directly inside of the standard frame should inherit from this class.'
};

sub _generate_content {
    my $self = shift;

    my $xml_view = $self->_get_xml_view(@_);

    my $xml_doc = XML::LibXML->createDocument();
    $self->_xml_doc($xml_doc);

    #these append to the XML doc
    $self->_generate_header();
    $self->_generate_body();
    $self->_generate_footer();

    my $xsl_doc = $self->_generate_xsl_doc($xml_view);
    return $self->transform_xml($xml_doc, $xsl_doc);
}

#override this for a custom title for the page
sub _title {
    my $self = shift;

    return $self->subject->__display_name__;
}

sub _generate_header {
    my $self = shift;
    my $xml_doc = $self->_xml_doc;

    # the header line is the class followed by the id
    my $object = $xml_doc->createElement('object');
    $object->addChild( $xml_doc->createAttribute('type', 'Genome::View::Static::Html'));
    $xml_doc->setDocumentElement($object);

    my $display_name = $object->addChild( $xml_doc->createElement('display_name') );
    $display_name->addChild( $xml_doc->createTextNode($self->_title) );

    return 1;
}

sub _generate_body {
    my $self = shift;
    my $xml_doc = $self->_xml_doc;

    my $html_text = '<body_html>' . $self->_html_body . '</body_html>';

    my $parser = XML::LibXML->new;
    my $html_xml_doc = $parser->parse_string($html_text);
    my $html_root = $html_xml_doc->documentElement;
    $xml_doc->adoptNode( $html_root );
    $xml_doc->documentElement->addChild( $html_root);

    return 1;
}

sub _html_body {
    my $self = shift;

    return '<div>To create a static view, implement _html_body (and optionally _title)</div>';
}

sub _generate_footer {
    my $self = shift;

    return 1;
}

sub _get_xml_view {
    my $self = shift;
    my %params = @_;

    my $xml_view = UR::Object::View->create(
        subject_class_name => $self->subject_class_name,
        perspective => 'default', #the desired_perspective is provided in _generate_xsl_doc
        toolkit => 'xml',
        %params,
    );

    return $xml_view;
}

sub transform_xml {
    my ($self,$xml_doc,$style_doc) = @_;

    my $parser = XML::LibXML->new;
    my $xslt = XML::LibXSLT->new;

    # convert the xml
    my $stylesheet = $xslt->parse_stylesheet($style_doc);
    my $results = $stylesheet->transform($xml_doc);
    my $content = $stylesheet->output_string($results);

    return $content;
}

sub _resolve_xsl_template_files {
    my $self = shift;
    my ($xml_view, $output_format, $xsl_path, $perspective) = @_;

    my @files = $self->SUPER::_resolve_xsl_template_files(@_);

    my $static_template = '/' . $self->output_format . '/static/root.xsl';
    push @files, $static_template
        if -f "$xsl_path$static_template";

    return @files;
}

1;
