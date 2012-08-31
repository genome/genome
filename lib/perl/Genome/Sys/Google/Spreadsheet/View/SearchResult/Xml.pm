package Genome::Sys::Google::Spreadsheet::View::SearchResult::Xml;

use strict;
use warnings;

use Genome;

class Genome::Sys::Google::Spreadsheet::View::SearchResult::Xml {
    is => 'Genome::View::SearchResult::Xml',
};


#rest_variable
#perspective
#subject
#toolkit
#solr_doc
#xsl_root



sub _generate_content {

    my ($self) = @_;

    my $doc = $self->solr_doc();

    return $doc->to_xml();
}




1;


