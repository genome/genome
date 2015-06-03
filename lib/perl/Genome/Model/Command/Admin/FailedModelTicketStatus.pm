package Genome::Model::Command::Admin::FailedModelTicketStatus;

use strict;
use warnings;

use Genome;

use Error qw(:try);
require RT::Client::REST;
require RT::Client::REST::Ticket;
require WWW::Mechanize;

BEGIN {
    $ENV{UR_DBI_NO_COMMIT} = 1;
}

class Genome::Model::Command::Admin::FailedModelTicketStatus {
    is => 'Genome::Command::WithColor',
    doc => 'report status of models in tickets',
    has_optional_input => [
        tickets => {
            is => 'Integer',
            doc => 'Which RT(s) to look at.  Defaults to all "new" and "open" tickets.',
            is_many => 1,
            shell_args_position => 1,
        },
        owner => {
            is => 'Text',
            doc => 'Limit any tickets to those owned by this user',
        },
        summary_only => {
            is => 'Boolean',
            doc => 'Only print status summary for each ticket',
            default => 0,
        },
    ],
};

sub help_detail {
    return <<HELP;
This command reports the status of models whose IDs are found in ticket(s).
HELP
}

sub _is_hidden_in_docs { return Genome::Sys->current_user_is_admin; }

sub execute {
    my $self = shift;

    # Connect
    my $rt = Genome::Model::Command::Admin::FailedModelTickets->_login_sso();

    my @ticket_ids = $self->tickets;
    unless(@ticket_ids) {
        @ticket_ids = Genome::Model::Command::Admin::FailedModelTickets::_find_open_tickets($self, $rt);
    }

    for my $ticket_id ( @ticket_ids ) {
        my $ticket = eval {
            RT::Client::REST::Ticket->new(
                rt => $rt,
                id => $ticket_id,
            )->retrieve;
        };
        unless ($ticket) {
            $self->error_message("Problem retrieving data for ticket $ticket_id: $@");
            next;
        }

        if(defined $self->owner and $ticket->owner ne $self->owner) {
            next;
        }

        my @ids;
        my $transactions = $ticket->transactions;
        my $transaction_iterator = $transactions->get_iterator;
        while ( my $transaction = &$transaction_iterator ) {
            my $content = $transaction->content;
            next unless $content;
            for my $line (split /\n/, $content) {
                my ($id) = $line =~ /(^[[:xdigit:]]{32})/;
                push @ids, $id if $id;
            }
        }

        my $owner = $ticket->owner;
        my $color = do {
            if    ($owner eq 'Nobody')              { 'yellow'  }
            elsif ($owner eq Genome::Sys->username) { 'green'   }
            else                                    { 'cyan'    }
        };

        $self->status_message(
            $self->_color($self->_color('Ticket %s', 'bold'), 'white') .' (' . $self->_color('%s',$color) . '):',
            $ticket_id, $ticket->owner
        );
        if (@ids) {
            my @models = Genome::Model->get(\@ids);
            if (@models) {
                Genome::Model::Command::Status->execute(models => \@models, summary_only => $self->summary_only);
            } else {
                $self->status_message('No models found for IDs.');
            }
        } else {
            $self->status_message('No model IDs found.');
        }
    }

    return 1;
}

1;
