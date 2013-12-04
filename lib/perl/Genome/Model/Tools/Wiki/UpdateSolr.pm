
package Genome::Model::Tools::Wiki::UpdateSolr;

use Genome;

use LWP::Simple;
use XML::Simple;


class Genome::Model::Tools::Wiki::UpdateSolr {
    is => ['Genome::Model::Tools::Wiki'],
    has => {
        days_ago => { 
            is => 'Text', 
            doc => 'How many days worth of wiki updates to send to solr', 
            default => '1',
            is_optional => 1 },
    },
    doc => 'Gets DAYS_AGO days worth of changes from wiki, submits to solr search engine for indexing',
};


sub help_synopsis {

    return "gmt wiki update-solr [ --days-ago N ]\n";
}


sub execute {

    my ($self) = @_;

    my $resource_lock = $ENV{GENOME_LOCK_DIR} . '/gcsearch/wiki_loader';

    if (Genome::Config->dev_mode()) {
        $resource_lock .= '_dev';
    }

    my $lock = Genome::Sys->lock_resource(resource_lock => $resource_lock, max_try => 1);
    die 'someone else has the wiki_loader lock' if !$lock;

    my $cache = Genome::Memcache->server() || die 'cant get memcache server';
    my $now = UR::Context->current->now();
    my $timeout = 60 * 60 * 24; # this is just storing which changes we've notified solr about
    print "inserting wiki pages to solr index at $now\n";

    # get/parse recent changes from wiki rss feed
    my $url = $self->url();
    my $raw_xml = LWP::Simple::get($url) || die "failed to get url: $url\n$!";
    my $parsed_xml = XML::Simple::XMLin($raw_xml);

    # title, link, description, pubDate, comments, dc:creator
    for my $item (@{ $parsed_xml->{'channel'}->{'item'} }) {

        # key is title and date of change
        my $key = cache_key($item);
        my $title = $item->{'title'};

        if ( ! defined($cache->get($key)) ) {

#            my $doc = Genome::Wiki::Document->get( title => $item->{'title'} )
#                || die 'cant get doc for title: ' . $item->{'title'};

            # NOTE: if this is slow we could Genome::Search->add(@all_docs)
            # but for now adding one, setting cache, adding another, setting cache...
            
            # post item to solr
#            Genome::Search->add($doc) || die 'Error: failed to add doc with title ' . $doc->title();
            Genome::Search::Queue->create(
                subject_class => 'Genome::Wiki::Document',
                subject_id => $title,
                priority => 9,
            );

            # mark as done
#            $cache->set($key, $now, $timeout ) || die "couldnt set cache for key=$key value=$now timeout=$timeout";
            # old version of Cache::Memcached was returning undef despite success
            $cache->set($key, $now, $timeout );

            printf("queued\t%s\t%s\n", $title, $key);
        }
    }

    UR::Context->commit();
    Genome::Sys->unlock_resource(resource_lock=>$lock);
}



sub cache_key {

    my ($item) = @_;

    my $title = $item->{'title'} || die 'couldnt make key- no title';
#    my $date = $item->{'pubDate'} || die "couldnt make key- no pubDate for item $title";

    my $key = join('---', 'wiki', $title);
    $key =~ s/\s//g; # memcache doesnt like whitespace

    return $key;
}



sub url {

    my ($self) = @_;

    my $url = join( '',
        $ENV{GENOME_SYS_SERVICES_WIKI_URL} . 'index.php?title=Special:RecentChanges&feed=rss&days=',
        $self->days_ago() );

    return $url;
}


1;


