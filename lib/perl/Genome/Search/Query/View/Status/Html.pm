package Genome::Search::Query::View::Status::Html;

use strict;
use warnings;
use Genome;

use Path::Class;

class Genome::Search::Query::View::Status::Html {
    is           => 'UR::Object::View::Default::Html',
    has_constant => [ perspective => { value => 'status', }, ],
};

## this is a cheap hack because i need to tell the search engine itself
## to give html, not transform it
sub _generate_content {

    my $self = shift;

    my $view = Genome::Search::Query::View::Status::Xml->create( subject => $self->subject(), format => 'xml' );
    my $raw_xml = $view->_generate_content();

    my $parser = XML::LibXML->new;
    my $xslt = XML::LibXSLT->new;
    my $source = $parser->parse_string($raw_xml);

    my $genomepm = file($INC{'Genome.pm'});
    my $root = $genomepm->dir();
    my $now = localtime();

    my $username = $ENV{'REMOTE_USER'};
    my $xsl_template = <<STYLE;
<?xml version="1.0" encoding="ISO-8859-1"?>
<xsl:stylesheet version="1.0"
xmlns:xsl="http://www.w3.org/1999/XSL/Transform">

<xsl:variable name="displayName">Search</xsl:variable>
<xsl:variable name="objectId">wee_id</xsl:variable>
<xsl:variable name="objectClassName">Genome::Search::Query</xsl:variable>
<xsl:variable name="currentTime">now</xsl:variable>
<xsl:variable name="currentPerspective">status</xsl:variable>
<xsl:variable name="currentToolkit">html</xsl:variable>
<xsl:variable name="username">$username</xsl:variable>
<xsl:variable name="noSearchBox">true</xsl:variable>
<xsl:variable name="resources">/view/genome/resource.html</xsl:variable>
<xsl:variable name="rest">/view</xsl:variable>

<xsl:include href="$root/Genome/xsl/html/status/root.xsl"/>
<xsl:include href="$root/Genome/xsl/html/status/ur_object.xsl"/>
<xsl:include href="$root/Genome/xsl/html/status/genome_search_query.xsl"/>
<xsl:include href="$root/Genome/xsl/html/search-result/universal.xsl"/>
<xsl:include href="$root/Genome/xsl/html/common.xsl"/>

<!-- <xsl:include href="$root/Genome/xsl/html/search-result/root.xsl"/> -->


</xsl:stylesheet>
STYLE

    my $style_doc = $parser->parse_string($xsl_template);
    my $stylesheet = $xslt->parse_stylesheet($style_doc);
    my $transformed_xml = $stylesheet->transform($source);
    my $content = $stylesheet->output_string($transformed_xml);

    return $content;
}




