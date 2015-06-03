package Genome::Model::Command::Admin::FailedModelTicketBase;

use strict;
use warnings;

use Error qw(:try);
require RT::Client::REST;
require RT::Client::REST::Ticket;
require WWW::Mechanize;

class Genome::Model::Command::Admin::FailedModelTicketBase {
    is => 'Genome::Command::WithColor',
    is_abstract => 1,
    doc => 'Base class for admin commands working with RT',
};

sub _find_open_tickets {
    my ($self, $rt) = @_;

    # The call to $rt->search() below messed up the login credentials stored in the
    # $rt session, making the loop at the bottom that retrieves tickets fail.
    # Save a copy of the login credentials here so we can re-set them when it's
    # time to get the ticket details
    my $login_cookies = $rt->_cookie();

    # Retrieve tickets -
    $self->status_message('Looking for tickets...');
    my @ticket_ids;
    try {
        @ticket_ids = $rt->search(
            type => 'ticket',
            query => "Queue = 'apipe-support' AND ( Status = 'new' OR Status = 'open' )",

        );
    }
    catch Exception::Class::Base with {
        my $msg = shift;
        if ( $msg eq 'Internal Server Error' ) {
            die 'Incorrect username or password';
        }
        else {
            die $msg->message;
        }
    };
    $self->status_message($self->_color('Tickets (new or open): ', 'bold').scalar(@ticket_ids));

    # re-set the login cookies that we saved away eariler
    $rt->_ua->cookie_jar($login_cookies);

    return @ticket_ids;
}

sub _ticket_for_id {
    my ($self, $rt, $ticket_id) = @_;

    my $ticket = eval {
        RT::Client::REST::Ticket->new(
            rt => $rt,
            id => $ticket_id,
        )->retrieve;
    };
    unless ($ticket) {
        $self->error_message("Problem retrieving data for ticket $ticket_id: $@");
    }

    return $ticket;
}

sub _server {
    return 'https://rt.gsc.wustl.edu/';
}

sub _login_sso {
    my $self = shift;

    my $mech = WWW::Mechanize->new(
        after =>  1,
        timeout => 10,
        agent =>  'WWW-Mechanize',
    );
    $mech->get( $self->_server() );

    my $uri = $mech->uri;
    my $host = $uri->host;
    if ($host ne 'sso.gsc.wustl.edu') {
        return;
    }

    $mech->submit_form (
        form_number =>  1,
        fields =>  {
            j_username => 'limsrt',
            j_password => 'Koh3gaed',
        },
    );
    $mech->submit();

    return RT::Client::REST->new(server => _server(), _cookie =>  $mech->{cookie_jar});
}

1;
