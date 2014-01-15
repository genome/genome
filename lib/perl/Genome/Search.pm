package Genome::Search;

use strict;
use warnings;

#Don't "use Genome;" here or we introduce a circular dependency.
use UR;
use MIME::Base64;

# JTAL: solr-dev is going to be prod, because old code will still point to solr

class Genome::Search {
    is => 'UR::Singleton',
    doc => 'This module contains methods for adding and updating objects in the Solr search engine.',
    has => [
        environment => {
            is => 'Text',
            calculate => q|
                return 'prod' if exists $ENV{GENOME_DEV_MODE} and $ENV{GENOME_DEV_MODE} == 1;
                return 'dev';
            |,
        },
        solr_server => {
            is => 'Text',
            default_value => $ENV{GENOME_SYS_SERVICES_SOLR},
            doc => 'Location of the Solr server',
        },
        _solr_server => {
            is => 'WebService::Solr',
            is_constant => 1,
            calculate_from => 'solr_server',
            calculate => q| require WebService::Solr; return WebService::Solr->new($solr_server); |,
        },
        cache_timeout => {
            is => 'Integer',
            default_value => 0,
            doc => 'Number of seconds for a document to persist in memcached.  Set to 0 for forever. [Note: If > 30 days, memcached instead uses the value as the timestamp at which the information should be expired.'
        },
        refresh_cache_on_add => {
            is => 'Boolean',
            default_value => 1,
            doc => 'If set, will cache the search result HTML when adding the item to the index.  If false, will clear the cache for the matching key but not update it.',
        },
        interpreter => {
            is => 'Text',
            calculate => q| return $^X; |,
        },
        genome_path  => {
            is => 'Text',
            calculate => q| return $INC{'Genome.pm'}; |,
        },
        ur_path => {
            is => 'Text',
            calculate => q| return $INC{'UR.pm'}; |,
        },
        workflow_path => {
            is => 'Text',
            calculate => q| return $INC{'Workflow.pm'}; |,
        }
    ],
};


#What classes are searchable is actually determined automatically by the existence of the relevant views.
#This just lists the order by which the results are typically sorted.
sub searchable_classes {

    my ($class) = @_;

    # order of this array determines sort order of search results

    my @ordered_searchable_classes = qw(
        Genome::Individual
        Genome::Taxon
        Genome::Library
        Genome::PopulationGroup
        Genome::Sample
        Genome::ModelGroup
        Genome::Model
        Genome::ProcessingProfile
        Genome::InstrumentData::FlowCell
        Genome::WorkOrder
        Genome::Site::TGI::Project
        Genome::Sys::Email
        Genome::DruggableGene::DrugGeneInteractionReport
        Genome::DruggableGene::DrugNameReport
        Genome::DruggableGene::GeneNameReport
        Genome::InstrumentData::Imported
        Genome::Sys::User
        Genome::Wiki::Document
    );

    return @ordered_searchable_classes;
}

sub environment {
    my $proto = shift;
    my $self = $proto->_singleton_object;

    if (@_ > 0) {
        $self->_solr_server(undef);
    }
    return $self->__environment(@_);
}

###  Index accessors  ###

sub search {
    my $class = shift;
    my $query = shift;
    my $webservice_solr_options = shift;

    my $self = $class->_singleton_object;
    my $response = $self->_solr_server->search($query, $webservice_solr_options);

    #TODO Better error handling--WebService::Solr doesn't handle error responses gracefully.
    return $response;
}

sub is_indexable {
    my $class = shift;
    my $object = shift;

    return $class->_resolve_solr_xml_view($object);
}

sub _resolve_solr_xml_view {
    my $class = shift;
    my $object = shift;

    my $subject_class_name = ref($object) || $object;

    return if (!$subject_class_name->isa('UR::Object'));
    return if $subject_class_name->isa('UR::Object::Ghost');
    return eval {
        UR::Object::View->_resolve_view_class_for_params(
            subject_class_name => $subject_class_name,
            toolkit => 'xml',
            perspective => 'solr'
        );
    };
}

sub _create_solr_xml_view_for_subject_class {
    my $class = shift;
    my $subject_class = shift;

    my $view_class = $class->_resolve_solr_xml_view($subject_class);
    unless ($view_class) {
        Carp::confess('To make an object searchable create an appropriate ::View::Solr::Xml that inherits from Genome::View::Solr::Xml.');
    }

    my $view = $view_class->create(subject_class_name => $subject_class, perspective => 'solr', toolkit => 'xml');
    return $view;
}

sub _resolve_result_xml_view {
    my $class = shift;
    my $subject_class_name = shift;

    return if (!$subject_class_name->isa('UR::Object'));
    return if $subject_class_name->isa('UR::Object::Ghost');
    return UR::Object::View->_resolve_view_class_for_params(
        subject_class_name => $subject_class_name,
        toolkit => 'xml',
        perspective => 'search-result'
    );
}


###  Index mutators  ###

sub add {
    my $class = shift;
    my @objects = @_;

    my $self = $class->_singleton_object;
    my @docs = $class->generate_document(@objects);

    unless($self->_solr_server->add(\@docs)) {
        $self->error_message('Failed to send ' . (scalar @docs) . ' document(s) to solr-dev.');
        return;
    }

#    my $solr_dev = WebService::Solr->new($self->_dev_solr_server_location());
#    unless($solr_dev->add(\@docs)) {
#        $self->error_message('Failed to send ' . (scalar @docs) . ' document(s) to other solr instance (solr)');
#        return;
#    }


    #$self->status_message('Sent ' . (scalar @docs) . ' document(s) to Solr.');
    return 1;
}

sub update {
    my $class = shift;

    #In solr, updating a record is the same as creating it--if the ID matches it overwrites
    return $class->add(@_);
}

sub delete {
    my $class = shift;
    my @objects = grep { exists $_->{db_committed} } @_;

    my $self = $class->_singleton_object;

    my @ids = map { join('---', $_->class, $_->id()) } @objects;

    my $error_count = $class->_delete_by_id(@ids);
    my $deleted_count = scalar(@ids) - $error_count;

    if($error_count) {
        $self->error_message('Failed to remove ' . $error_count . ' document(s) from Solr.');
        return;
    }

    return $deleted_count || 1;
}

sub delete_by_class_and_id {
    my $class = shift;
    my $subject_class = shift;
    my $subject_id = shift;

    my $self = $class->_singleton_object;

    my $solr = $self->_solr_server;

    my $memcached = Genome::Memcache->server;

    my $solr_index_id = join('---', $subject_class, $subject_id);
    my $cache_id = join('genome_search:', $solr_index_id);

    if ($solr->delete_by_id($solr_index_id)) {
        $memcached->delete($cache_id);
    }
    else {
        $self->error_message("Failed to remove document from search (Class: $subject_class, ID: $subject_id).");
        return;
    }

    return 1;
}

sub _delete_by_id {
    my ($class, @ids) = @_;

    my $self = $class->_singleton_object;

    my $solr = $self->_solr_server();

    my $memcached = Genome::Memcache->server;

    my $error_count = 0;
    for my $id (@ids) {
        if($solr->delete_by_id($id)) {
            my $cache_key = "genome_search:$id";
            $memcached->delete($cache_key);
        } else {
            $error_count++;
        }
    }

    return $error_count;
}

sub delete_by_doc {
    my $class = shift;
    my @docs = @_;

    unless (@docs) {
        return;
    }

    my $self = $class->_singleton_object;
    my $solr = $self->_solr_server;
    my $memcached = Genome::Memcache->server;

    my (@cache_ids, @solr_ids);
    for my $doc (@docs) {
        my $solr_id = $doc->value_for('id');
        unless ($solr_id) {
            die $self->error_message('No value for \'id\':' . Dumper($doc));
        }
        push @solr_ids, $solr_id;
        push @cache_ids, "genome_search:" . $solr_ids[-1];
    }

    my $response = $solr->_send_update('<delete><id>' . join('</id><id>', @solr_ids) . '</id></delete>');
    if ($response->ok) {
        for my $cache_id (@cache_ids) {
            $memcached->delete($cache_id);
        }
    }

    return $response;
}

sub clear {
    my $class = shift;

    my $self = $class->_singleton_object;

    return 1 if UR::DBI->no_commit; #Prevent automated index manipulation when changes certainly won't be committed

    my $solr = $self->_solr_server;

    $solr->delete_by_query('*:*') || return; #Optimized by solr for fast index clearing
    $solr->optimize() || return; #Prevent former entries from influencing future index

    #$self->status_message('Solr index cleared.');

    #NOTE: The memcached information is not cleared at this point.
    #However, anything added to search will trigger a cache update.

    return 1;
}

sub cache_key_for_doc {
    my $class = shift;
    my $doc = shift;

    return 'genome_search:' . $doc->value_for('id');
}


### Other ###

sub get_subject_from_doc {
    my $class = shift;
    my $solr_doc = shift || die;

    my $subject_class_name = $solr_doc->value_for('class');
    return unless $subject_class_name->can('get');

    my $object_id = $solr_doc->value_for('object_id');
    if ($object_id) {
        if ($subject_class_name->isa('UR::Object::Set') && $object_id =~ /=$/) {
            # Sets have invalid XML chars in their IDs so we encode them in Base64.
            # Encoding is done in Genome::View::Solr::Xml::_generate_object_id_field_data so keep symmetry there.
            # TODO Encoding/decoding should probably be handled by the object itself.
            $object_id = decode_base64($object_id);
        }
        return $subject_class_name->get($object_id);
    }

    my $solr_id = $solr_doc->value_for('id');
    if ($solr_id) {
        my ($derived_object_id) = $solr_id =~ /.*?(\d+)$/;
        if ($derived_object_id) {
            return $subject_class_name->get($derived_object_id);
        }
    }

    return;
}

###  XML Generation for results  ###

sub generate_pager_xml {
    my $class = shift;
    my $pager = shift;
    my $xml_doc = shift || XML::LibXML->createDocument();

    my $page_info_node = $xml_doc->createElement('page-info');

    $page_info_node->addChild( $xml_doc->createAttribute('previous-page', $pager->previous_page) )
        if $pager->previous_page;
    $page_info_node->addChild( $xml_doc->createAttribute('current-page', $pager->current_page) );
    $page_info_node->addChild( $xml_doc->createAttribute('next-page', $pager->next_page) )
        if $pager->next_page;
    $page_info_node->addChild( $xml_doc->createAttribute('last-page', $pager->last_page) );

    return $page_info_node;
}

sub generate_result_xml {
    my $class = shift;
    my $doc_or_docs = shift;
    my $xml_doc = shift || XML::LibXML->createDocument();
    my $format = shift || 'xml';
    my $override_cache = shift || 0;

    my @docs;
    if(ref $doc_or_docs eq 'ARRAY') {
        @docs = @$doc_or_docs;
    } else {
        @docs = $doc_or_docs;
    }

    my %views;
    my @result_nodes;

    for my $doc (@docs) {
        my $object_id = $doc->value_for('object_id');
        my $object_class = $doc->value_for('class');

        unless($object_id) {
            #Older snapshots will create duplicate entries in the index; let's not show them.
            $class->delete_by_doc($doc);
            next;
        }

        if(!$override_cache and $format eq 'html') {
            if(my $result_node = $class->_get_cached_result($doc, $xml_doc)) {
                push @result_nodes, $result_node;
                next;
            }
        }


        require Genome; #It is only at this point that we actually need to load other objects

        $object_class->can('isa'); #Force class autoloading
        my $object = $object_class->get($object_id);

        unless($object and ($object_class eq $object->class)) {
            #Either
            # (1) the entity in the index no longer exists
            # (2) the entity has a different class than was indexed(!) so index is wrong
            #so remove the object from the index
            #(in the case of #2 the correct class will be indexed by the cron later)
            $class->delete_by_doc($doc);
            next;
        }

        my $view;
# NOTE: turning off this optimization; it reuses the first object, which doesnt work
#  out so well when your view is doing $self->property (self is the original instance)


        if($views{$object->class}) {
            $view = $views{$object->class};
            $view->subject($object);
            $view->solr_doc($doc);
            $view->_update_view_from_subject();
        } else {
            if(my $view_class = $class->_resolve_result_xml_view($object_class)) {
                my %view_args = (
                    perspective => 'search-result',
                    toolkit => $format,
                    solr_doc => $doc,
                    subject => $object,
                    rest_variable => '/view',
                );

                if($format eq 'xsl' or $format eq 'html') {
                    $view_args{xsl_root} = Genome->base_dir . '/xsl';
                }

                $view = $view_class->create(%view_args);
                $views{$object_class} = $view;
            } else {
                $class->error_message('No suitable search result view found for ' . $object_class . '.');
                next;
            }
        }


        my $object_content = $view->content;

        if($format eq 'xsl' or $format eq 'html') {
            $object_content =~ s/^<\?.*?\?>//;
        }

        my $result_node = $xml_doc->createElement('result');

        if($format eq 'xml') {
            my $lib_xml = XML::LibXML->new();
            my $content = $lib_xml->parse_string($object_content);

            $result_node->addChild($content->childNodes);
        } else {
            $result_node->addChild($xml_doc->createTextNode($object_content));
        }

        push @result_nodes, $result_node;

        if($format eq 'html') {
            $class->_cache_result($doc, $result_node);
        }
    }

    return @result_nodes;
}

sub _cache_result {
    my $class = shift;
    my $doc = shift;
    my $result_node = shift;

    my $self = $class->_singleton_object;

    my $html_to_cache = $result_node->childNodes->string_value;

    my $memcached = Genome::Memcache->server;

    return $memcached->set($self->cache_key_for_doc($doc), $html_to_cache, $self->cache_timeout);
}

sub _delete_cached_result {
    my $class = shift;
    my $doc = shift;

    my $memcached = Genome::Memcache->server;

    $memcached->delete($class->cache_key_for_doc($doc));
}

sub _get_cached_result {
    my $class = shift;
    my $doc = shift;
    my $xml_doc = shift;

    my $memcached = Genome::Memcache->server;

    my $cache_key = $class->cache_key_for_doc($doc);
    my $html_snippet = $memcached->get($cache_key);

    #Cache miss
    return unless $html_snippet;

    my $result_node = $xml_doc->createElement('result');
    $result_node->addChild($xml_doc->createTextNode($html_snippet));

    return $result_node;
}

###  Search "document" creation/delegation  ###

my %views;
sub generate_document {
    my $class = shift;
    my @objects = @_;

    my @docs;
    for my $o (@objects) {
        # Building new instances of View classes is slow as the system has to resolve a large set of information
        # According to NYTProf, recycling the view reduced the time spent here from 457s to 45.5s on a set of 1000 models
        my $view = $views{$o->class} || $class->_create_solr_xml_view_for_subject_class($o->class);
        $view->subject($o);
        $view->_update_view_from_subject();
        push @docs, $view->content_doc;
    }
    return @docs;
}

###  Callback for automatically updating index  ###

our $LOADED_MODULES = 0;
sub _index_queue_callback {
    my ($class, $object, $aspect) = @_;

    return unless $object;

    unless ($LOADED_MODULES++) {
        require WebService::Solr;
        require MRO::Compat;
    }

    my $index_queue;
    my $meta = $object->__meta__;
    my @trigger_properties = ('create', 'delete', $meta->all_property_names);
    if (grep { $aspect eq $_ } @trigger_properties) {
        $index_queue = Genome::Search->queue_for_update($object);
    }

    return $index_queue;
}

sub queue_for_update {
    my $class = shift;
    my @objects = @_;

    return unless @objects;

    my @index_queues;
    for my $object (@objects) {
        my %create_params = (
            subject_id => $object->id,
            subject_class => $object->class,
        );
        if ($object->class->can('search_index_queue_priority')) {
            $create_params{priority} = $object->class->search_index_queue_priority;
        }
        my $index_queue = Genome::Search::Queue->create(%create_params);
        push @index_queues, $index_queue;
    }

    if(@objects == 1) {
        return $index_queues[0];
    } else {
        return @index_queues;
    }
}

my $observer;
# register_callbacks and unregister_callbacks should be called from Genome.pm,
# so typically it won't need to be called elsewhere.
sub register_callbacks {
    my $class = shift;
    my $searchable_class = shift;

    $observer = $searchable_class->add_observer(
        callback => sub { $class->_index_queue_callback(@_); },
    );
}

sub unregister_callbacks {
    $observer->delete unless $observer->isa("UR::DeletedRef");
}

sub get_indexed_document_count {
    my $class = shift;
    my $self = $class->_singleton_object();

    my $response = $class->search(
        '*:*',
        {
            rows => 0,
            hl => 'false',
            defType => 'lucene'
        }
    );

    if ($response->ok) {
        return $response->content->{response}{numFound};
    } else {
        die($self->error_message($response->status_message));
    }
}

sub get_indexed_class_counts {
    my $class = shift;
    my $self = $class->_singleton_object();

    my $response = $class->search(
        '*:*',
        {
            rows => 0,
            hl => 'false',
            defType => 'lucene',
            'facet.field' => 'class'
        }
    );

    if ($response->ok) {
        return @{$response->content->{facet_counts}{facet_fields}{class}};
    } else {
        die($self->error_message($response->status_message));
    }
}

#OK!
1;

=pod

=head1 NAME

Genome::Search

=head1 SYNOPSIS

  Genome::Search->add(@objects);
  Genome::Search->delete(@objects);
  Genome::Search->is_indexable($object);

  Genome::Search->search($query, $options);

=head1 DESCRIPTION

This class adds, updates, and deletes entries for objects from the Solr index.

=head1 METHODS

=over 4

=item search

Query the Solr index.

=item is_indexable

Determine if an appropriate subclass exists to add the object to the Solr index.

=item add

Adds one or more objects to the Solr index.

=item update

An alias for add. Solr's add method automatically overwrites any existing entry for the object.

=item delete

Removes one or more objects from the Solr index.

=item clear

Removes all objects from the Solr index.

=back

=head1 DETAILS

To make a new type of object indexable by Solr, create a subclass of Genome::View::Solr::Xml for
that object.  For example, to make a Genome::Individual indexable, the class
Genome::Individual::View::Solr::Xml was created.  The class should define the
"type" attribute with an alphanumeric name representing the class in the index
as well as specify what aspects to be included in what search fields.

Additionally, for the display of results, appropriate ::View::SearchResult::Xml and
::View::SearchResult::Html classes should be created for each class where a
::View::Solr::Xml has been defined.

=head1 SEE ALSO

Genome::Individual::View::Solr::Xml
Genome::Individual::View::SearchResult::Xml
Genome::Individual::View::SearchResult::Html

=cut
