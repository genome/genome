package Genome::Search::Query::View::Status::Xml;

use strict;
use warnings;
use Genome;

use XML::Simple;

class Genome::Search::Query::View::Status::Xml {
    is           => 'UR::Object::View::Default::Xml',
    has_constant => [ perspective => { value => 'status', }, ],
};

my $RESULTS_PER_PAGE = 20;

sub _generate_content {
    my $self    = shift;
    my $subject = $self->subject;
    my $format = exists $subject->{format} ? $subject->{format} : 'xml';

    my $query = $subject->query;
    my $page  = $subject->page;

    my $doc          = XML::LibXML->createDocument();
    my $results_node = $doc->createElement('solr-results');

    my $params = {   
        rows         => $RESULTS_PER_PAGE,
        start        => $RESULTS_PER_PAGE * ( $page - 1 )
    };

    my $fq = $subject->fq;    # facet query -- see dismax docs
    if ($fq) {
        # during drill down, still show all facet counts
        $params->{'facet.field'} = '{!ex=tt}type';
        $params->{'fq'} = '{!tag=tt}' . $fq;
   
        my ($facet_name) = $fq =~ /type:(.*)/;
        $facet_name =~ s/"//g;
        $results_node->addChild( $doc->createAttribute( "facet-name", $facet_name ));
    }

    $params->{'hl'} = 'true';
#    $params->{'hl.fl'} = 'title,content';

    $params->{'qs'} = 1;

    my $solrQuery = $query;
    my $response = Genome::Search->search(
        $solrQuery,
        $params
    );

    my @p;
    for my $k (keys %$params) {
        push @p, join('', $k, '=', $params->{$k});
    }

    my $param_str = join(',', @p);

    my $time = $UR::Context::current->now;
    $results_node->addChild( $doc->createAttribute( "generated-at", $time ) );
    $results_node->addChild( $doc->createAttribute( "input-name",   "query" ) );
    $results_node->addChild( $doc->createAttribute( "query",        $query ) );
    $results_node->addChild( $doc->createAttribute( "params",        $param_str) );
    $results_node->addChild( $doc->createAttribute( "num-found", $response->content->{'response'}->{'numFound'} ));
#    $results_node->addChild( $doc->createAttribute( "auth-user", $ENV{'REMOTE_USER'} || Genome::Sys->username() ));

#   HIGHLIGHTING XML

    my $highlights_node = $doc->createElement('highlights');

    my $highlights_raw = $response->content->{'highlighting'};
    for my $found_id (keys %$highlights_raw) {
        my $h = $highlights_raw->{$found_id};
        my $field_str;
        for my $field (keys %$h) {
            my @field_highlights = @{ $h->{$field} };
            $field_str = "$field: " . join(',', @field_highlights) . '</br>';
        }
        my $highlight_item = $doc->createElement('highlight_item');
        $highlight_item->addChild($doc->createAttribute('id',$found_id));
        $highlight_item->addChild($doc->createAttribute('description',$field_str));
        $highlights_node->addChild($highlight_item);
    }

    $results_node->addChild($highlights_node);


#   END OF HIGHLIGHTING XML


#   FACET XML
#    facet_dates
#    facet_fields
#    facet_queries

    my $f = $response->facet_counts();
        
    if ($f) {
    
        my $facets_node = $doc->createElement('facets');

        my @raw_fields = @{ $f->{'facet_fields'}->{'type'} };
        my $facets = {};

        for (my $i = 0; $i < scalar(@raw_fields); $i+=2) {
            $facets->{ $raw_fields[$i] } = $raw_fields[$i + 1];
        }

        my $facet_total = 0;

        for my $field_name (sort keys %$facets) {
            my $count = $facets->{$field_name};
            next if ! $count;
            my $field = $doc->createElement('field');
            $field->addChild( $doc->createAttribute('name', $field_name) );
            $field->addChild( $doc->createAttribute('label', make_label($field_name)) );
            $field->addChild( $doc->createAttribute('icon-prefix', icon_prefix($field_name)) );
            $field->addChild( $doc->createAttribute('count', $count) );
            $facets_node->addChild($field);

            $facet_total += $count;
        }
        
        $results_node->addChild($doc->createAttribute('facet-total', $facet_total));
        $results_node->addChild($facets_node);
    }


#   END OF FACET XML GENERATION

    # create query-no-types attribute
    my @params = split /\s+/, $query;
    for ( my $i = $#params ; $i >= 0 ; --$i ) {
        splice @params, $i, 1
          if $params[$i] =~ /(type\:(\S+))/i;
    }
    my $query_no_types = join " ", @params;

    $results_node->addChild(
        $doc->createAttribute( "query-no-types", $query_no_types ) );
    $results_node->addChild(
        Genome::Search->generate_pager_xml( $response->pager, $doc ) );

    # if using sort_solr_docs(), each page is sorted but overall
    # the results are not sorted
#    my @ordered_docs = sort_solr_docs( $response->docs );
    my @ordered_docs = $response->docs();

#    my @result_nodes =
#      Genome::Search->generate_result_xml( \@ordered_docs, $doc, $format );
    my @docs = @ordered_docs;

#            my $raw_xml = XMLout($doc, 'root_name' => 'result');
#            $result_node = XML::LibXML->load_xml(string => $raw_xml);

    for my $doc (@docs) {
        my $raw_xml = $doc->to_xml();
        my $result_node = XML::LibXML->load_xml(string => $raw_xml);
        $results_node->addChild($result_node->documentElement->cloneNode(1));
    }

    $doc->setDocumentElement($results_node);

    $doc->toString(1);
}


sub icon_prefix {

    # "type" field in solr index needs to be cleaned up;
    # for now, we have this hash to map to icon url
    
    my ($type) = @_;

    my $icon_prefix = {
        'model'                  => 'genome_model',
        'model - lane_qc'        => 'genome_model',
        'model - alignment'      => 'genome_model',
        'model - microarray'     => 'genome_model',
        'model - somatic'        => 'genome_model',
        'model - other'          => 'genome_model',
        'model - rna'            => 'genome_model',
        'work-order'             => 'genome_workorder',
        'project'                => 'genome_project',
        'mail'                   => 'genome_sys_email',
        'wiki-page'              => 'genome_wiki_document',
        'modelgroup'             => 'genome_modelgroup',
        'population_group'       => 'genome_populationgroup',
        'disk_group'             => 'genome_disk_group',
        'disk_volume'            => 'genome_disk_volume',
        'illumina_run'           => 'genome_instrumentdata_flowcell',
        'individual'             => 'genome_individual',
        'library'                => 'genome_library',
        'processing_profile'     => 'genome_processingprofile',
        'sample'                 => 'genome_sample',
        'taxon'                  => 'genome_taxon',
        'user'                   => 'genome_sys_user',
        'solexa_instrument_data' => 'genome_instrumentdata',
    };

    return $icon_prefix->{$type};

}

sub make_label {

    my ($text) = @_;

    $text =~ s/[_]/ /g;

    if ($text =~ /^Model\s\w+/) {
        $text =~ s/^Model/Model -/;
    }

    my @words;
    for my $w (split(/\s+/,$text)) {
        if($w eq 'qc') {
            push @words, 'QC';
        } else {
            push @words, ucfirst($w);
        }
    }

    return join(' ',@words);
}

sub sort_solr_docs {
    my @docs = @_;

    my @ordered_doc_classes = Genome::Search->searchable_classes();

    my %ordered_docs;
    my @everything_else_docs;

    my %docs_by_class;

    for my $solr_doc (@docs) {
        my $this_doc_class = $solr_doc->value_for('class');
        my ($matched_class) =
          grep { $this_doc_class =~ m/$_/ } @ordered_doc_classes;
        if ($matched_class) {
            push @{ $docs_by_class{$matched_class} }, $solr_doc;
        } else {
            push @everything_else_docs, $solr_doc;
        }
    }

    my @doc_classes = keys %docs_by_class;
    my @ordered_docs;

    for my $ordered_class (@ordered_doc_classes) {
        push @ordered_docs, @{ $docs_by_class{$ordered_class} }
          if ( exists $docs_by_class{$ordered_class} );
    }
    push @ordered_docs, @everything_else_docs;

    return @ordered_docs;
}

