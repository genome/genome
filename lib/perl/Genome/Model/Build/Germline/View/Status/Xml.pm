#:boberkfe this looks like a good place to use memcache to cache up some build status.
#:boberkfe when build events update, stuff their status into memcache.  gathering info otherwise
#:boberkfe can get reaaaaal slow.

package Genome::Model::Build::Germline::View::Status::Xml;

use strict;
use warnings;
use Genome;
use Data::Dumper;
use XML::LibXML;
use XML::LibXSLT;

class Genome::Model::Build::Germline::View::Status::Xml {
    is => 'Genome::Model::Build::View::Status::Xml',
};

sub get_build_node {
    my $self = shift;

    my $node = $self->SUPER::get_build_node(@_);
    my $doc = $self->_doc;
    my $build = $self->subject;

    my $base_build = $build->source_build;
    my $accumulated_alns = $base_build->accumulated_alignments_directory;

    $node->addChild($doc->createAttribute('base-alignment-path', $accumulated_alns));

    return $node;
}

1;
