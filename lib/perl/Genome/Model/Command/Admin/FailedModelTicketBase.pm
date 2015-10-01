package Genome::Model::Command::Admin::FailedModelTicketBase;

use strict;
use warnings;

use Error qw(:try);
require RT::Client::REST;
require RT::Client::REST::Ticket;
require WWW::Mechanize;

class Genome::Model::Command::Admin::FailedModelTicketBase {
    is => 'Command::V2',
    roles => 'Genome::Role::CommandWithColor',
    is_abstract => 1,
    doc => 'Base class for admin commands working with RT',
};

sub _find_open_tickets {
    my ($self, $rt) = @_;

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

sub _login_sso {
    my $self = shift;

    my $mech = WWW::Mechanize->new(
        after =>  1,
        timeout => 10,
        agent =>  'WWW-Mechanize',
    );
    $mech->get( Genome::Config::get('rt_url') );

    my $uri = $mech->uri;
    my $host = $uri->host;
    if ($host ne 'sso.gsc.wustl.edu') {
        return;
    }

    $mech->submit_form (
        form_number =>  1,
        fields =>  {
            j_username => Genome::Config::get('rt_login'),
            j_password => Genome::Config::get('rt_auth'),
        },
    );
    $mech->submit();

    my $rt = RT::Client::REST->new(server => Genome::Config::get('rt_url'), _cookie =>  $mech->{cookie_jar});
    #propagate the cookie jar to the UA since we're not calling $rt->login
    $rt->_ua->cookie_jar( $mech->{cookie_jar} );

    return $rt;
}

1;
