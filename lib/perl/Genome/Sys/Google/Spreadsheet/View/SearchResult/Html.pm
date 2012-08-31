package Genome::Sys::Google::Spreadsheet::View::SearchResult::Html;

use strict;
use warnings;

use Genome;

class Genome::Sys::Google::Spreadsheet::View::SearchResult::Html {
    is => 'Genome::View::SearchResult::Html',
};


sub _generate_content {

    my ($self) = @_;

    return $self->SUPER::_generate_content(solr_doc => $self->solr_doc());
}

1;
