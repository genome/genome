
package Genome::Wiki::Document;

use Genome;
use LWP::UserAgent;
use XML::Simple;
use HTML::Parser;

use Carp;

# why are some of these things calculated values?
# because default_value wasnt working as expected-
# the other option was to make methods below that calle
# __myProperties

class Genome::Wiki::Document {
    is  => ['UR::Value','Genome::Searchable'],
    id_by => [
        title => { is => 'Text' },
    ],
    doc => 'represents a page in the wiki',
    has_transient => {
        revision_id   => { is => 'Number' },
        content       => { is => 'Text' },
        timestamp     => { is => 'Text' },
        user          => { is => 'Text' },
    },
    has => {
        environment => {
             calculate => q{ Genome::Config::dev_mode() ? 'dev' : 'prod' },
        },
        wiki_server_url => {
            calculate => qq{ '$ENV{GENOME_SYS_SERVICES_WIKI_URL}' . 'api.php' },
        },
    },
};


# This class is a UR::Value- when you do a get, it pulls the page from the wiki and
# represents it in a Genome::Wiki::Document


__PACKAGE__->add_observer (
    aspect   => 'load',
    callback => sub {
        my ($self) = @_;
        $self->initialize();
    }
);


sub initialize {

    my ($self) = @_;

    my $ua        = LWP::UserAgent->new();
    my $query_url = $self->query_url();

    my $r         = $ua->get($query_url);

    if (! $r->is_success ) {
        confess('couldnt get wiki doc: ' . $self->query_url());
    }

    my $xml = $r->content;
    $self->absorb_xml($xml);

    return 1;
}

sub absorb_xml {

    my ($self, $xml) = @_;

    my $p = XMLin($xml) || confess('failed to parse xml doc');

    my $page = $p->{'query'}->{'pages'}->{'page'}->{'revisions'}->{'rev'};

    my $content = strip_html_from_string($page->{'content'});
    $self->content($content);
    $self->revision_id($page->{'revid'});
    $self->timestamp($page->{'timestamp'});
    $self->user($page->{'user'});

    return 1;
}

sub query_url {
    
    my ($self) = @_;

    my $server_url = $self->wiki_server_url();
    my $query = $self->build_query();

    my $params = join( '&', 'action=query', 
                            'format=xml', 
                            'prop=revisions',
                            'rvprop=content|ids|user|timestamp',
                            $query
    );

    my $query_url = join('?', $server_url, $params);

    return $query_url;
}

sub build_query {
    
    my ($self) = @_;
    my $q;

    if ($self->title()) {
        $q = join('=', 'titles', $self->title()); 
    } elsif ($self->revision_id()) {
        $q = join('=', 'revids', $self->revision_id());
    } else {
        confess('cant build query- require a title or revision id');
    }

    return $q;
}

sub strip_html_from_string {

    my ($s) = @_;

    my %inside;
    my $stripped_content;

    my $tag_function = sub {
        my ( $tag, $num ) = @_;
        $inside{$tag} += $num;
        $stripped_content .= " ";    # not for all tags
    };

    my $text_function = sub {
        return if $inside{script} || $inside{style};
        $stripped_content .= $_[0];
    };

    # TODO: not passing stuff in at function
    HTML::Parser->new(
        handlers => [
            start => [ $tag_function,  "tagname, '+1'" ],
            end   => [ $tag_function,  "tagname, '-1'" ],
            text  => [ $text_function, "dtext" ],
        ],
        marked_sections => 1,
        )->parse($s)
        || die 'couldnt parse for stripping out html from content';

    return $stripped_content;
}


1;



