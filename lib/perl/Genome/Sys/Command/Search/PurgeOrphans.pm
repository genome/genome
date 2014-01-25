package Genome::Sys::Command::Search::PurgeOrphans;

use strict;
use warnings;

use Genome;
use Genome::Utility::PluckColumn;
use WebService::Solr::Query;

class Genome::Sys::Command::Search::PurgeOrphans {
    is => 'Command::V2',
    doc => 'Remove documents from search index which no longer have corresponding objects in the database.',
};

sub help_detail {

return <<EOS
Scans the search engine index for documents which no longer have corresponding
objects in the database.

Theoretically this should never catch any to remove but in practice we have
seen it trigger in small amounts.  The current belief is that those are due to
somebody removing the object after the initial Solr query but before that
object is checked in this process. This also means that their orphan state is a
false positive.
EOS

}

sub execute {
    my $self = shift;

    my @class_names = keys %{{Genome::Search->get_indexed_class_counts()}};
    for my $class_name (@class_names) {

        my $pid = UR::Context::Process->fork();
        if (not defined $pid) {
            die "Failed to fork ($!).";
        }
        elsif ($pid == 0) {
            # CHILD
            my $start_time = time;
            my @docs = $self->get_documents_for_class($class_name);
            my @docs_to_purge = $self->find_orphaned_docs($class_name, @docs);
            $self->debug_message("Processed " . @docs . " $class_name documents in " . (time - $start_time) . " seconds and found " . @docs_to_purge . " to remove.");

            if (@docs_to_purge) {
                $self->debug_message("Example $class_name id being purged: %s", $docs_to_purge[0]->value_for('object_id'));
                $start_time = time;
                my $response = Genome::Search->delete_by_doc(@docs_to_purge);
                if ($response->ok) {
                    $self->status_message("Removed " . @docs_to_purge . " documents in " . (time - $start_time) . " seconds:");
                    $self->status_message("\t" . join("\n\t", map { $_->value_for('id') } @docs_to_purge) . "\n");
                } else {
                    $self->status_message("Failed to remove " . @docs_to_purge . " documents.");
                }
            }

            exit;
        }
        else {
            # PARENT
            waitpid($pid, 0);
        }
    }
    return 1;
}

sub get_documents_for_class {
    my $self = shift;
    my $class_name = shift;

    die("Must supply a class name!") unless $class_name;

    my $response = Genome::Search->search(
        'id:[* TO *]',
        {
            start => 0,
            rows => Genome::Search->get_indexed_document_count(),
            fl => 'id,object_id',
            hl => 'false',
            defType => 'lucene',
            fq => [
                WebService::Solr::Query->new( { class => $class_name } )
            ]
        },
    );
    unless($response->ok) {
        die($self->error_message("Search response failed, aborting."));
    }
    return $response->docs;
}

sub find_orphaned_docs {
    my $self = shift;
    my $class_name = shift;
    my @docs = @_;

    my %ids_in_database_lookup;
    my @docs_to_purge;

    my @id_properties = grep $_, $class_name->__meta__->get_all_id_column_names;

    if (@id_properties > 1) {
        $self->debug_message("Composite IDs are not currently supported! Class: $class_name");
        return ();
    }
    if (scalar(@id_properties) == 0) {
        $self->debug_message("Skipping $class_name - its not database backed!");
        return ();
    }

    my $ids_in_database = Genome::Utility::PluckColumn::pluck_column_from_class($class_name, column_name => $id_properties[0]);
    $ids_in_database_lookup{$_} = 1 for(@$ids_in_database);

    for my $current_doc (@docs) {
        my $obj_id = $current_doc->value_for('object_id');
        push @docs_to_purge, $current_doc unless $ids_in_database_lookup{$obj_id};
    }

    return @docs_to_purge;

}

1;
