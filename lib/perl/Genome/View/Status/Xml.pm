package Genome::View::Status::Xml;

use strict;
use warnings;
use Genome;

class Genome::View::Status::Xml {
    is => 'UR::Object::View::Default::Xml',
#    is_abstract => 1,
    has_constant => [
        perspective => {
            value => 'status',
        },
    ],
    doc => 'The base class for creating the XML document representing the full-page status of an object'
};

sub _generate_content {
    my $self = shift;

    # classes inherit from this, but i'm not sure thats a good idea
    # this should workaround any problems caused by it
    unless ($self->subject->isa('UR::Namespace')) {
        return $self->SUPER::_generate_content;
    }

    my $doc = XML::LibXML->createDocument();
    my $search_form_node = $doc->createElement("search-form");
    my $time = $UR::Context::current->now;

    $search_form_node->addChild( $doc->createAttribute("generated-at",$time) );

    $doc->setDocumentElement($search_form_node);

    # get list of events
    my $sql = "SELECT DISTINCT event_status FROM mg.genome_model_event ORDER BY event_status";

    my $dbh = Genome::DataSource::GMSchema->get_default_handle();
    my $sth = $dbh->prepare($sql) or die "Failed to connect to the database!";
    $sth->execute() or die "Error querying the database!";

    # create and populate event-statuses node
    my $events_node = $search_form_node->addChild( $doc->createElement("event-statuses") );
    while (my @events_data = $sth->fetchrow_array()) {
        if($events_data[0]) {
            my $event_node = $events_node->addChild( $doc->createElement("event-status") );
            $event_node->addChild( $doc->createTextNode($events_data[0]) );
        }
    }

    # get list of users
    my $sql2 = "SELECT DISTINCT user_name FROM mg.genome_model_event ORDER BY user_name";

    my $dbh2 = Genome::DataSource::GMSchema->get_default_handle();
    my $sth2 = $dbh2->prepare($sql2) or die "Failed to connect to the database!";
    $sth2->execute() or die "Error querying the database!";

    # create and populate users node
    my $users_node = $search_form_node->addChild( $doc->createElement("users") );
    while (my @users_data = $sth2->fetchrow_array()) {
        if ($users_data[0]) {
            my $user_node = $users_node->addChild( $doc->createElement("user") );
            $user_node->addChild( $doc->createTextNode($users_data[0]) );
        }
    }

    $doc->toString(1);
}

1;
