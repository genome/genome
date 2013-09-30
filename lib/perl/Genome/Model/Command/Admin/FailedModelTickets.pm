package Genome::Model::Command::Admin::FailedModelTickets;

use strict;
use warnings;

use Genome;

use Data::Dumper 'Dumper';
use Error qw(:try);
use File::Find 'find';
use File::Grep 'fgrep';
require IO::Prompt;
require RT::Client::REST;
require RT::Client::REST::Ticket;
require WWW::Mechanize;

class Genome::Model::Command::Admin::FailedModelTickets {
    is => 'Command::V2',
    doc => 'find failed cron models, check that they are in a ticket',
    has_input => [
        include_failed => {
            is => 'Boolean',
            default_value => 1,
            doc => 'Include builds with status Failed',
        },
        include_unstartable => {
            is => 'Boolean',
            default_value => 1,
            doc => 'Include builds with status Unstartable',
        },
        include_pending => {
            is => 'Boolean',
            default_value => 0,
            doc => 'Include builds whose model status is requested, scheduled, or running.',
        },
    ],
};

sub help_detail {
    return <<HELP;
This command collects cron models by failed or unstartable build events and scours tickets for them. If they are not found, the models are summaraized first by the error entry log and then by grepping the error log files. The summary is the printed to STDOUT.
HELP
}

sub execute {
    my $self = shift;

    my @events = $self->get_events();
    my %builds = $self->get_builds(@events);

    my %tickets = $self->get_tickets(%builds);
    for my $model_id (keys %tickets) {
        delete $builds{$model_id};
    }

    # Consolidate errors
    $self->status_message('Consolidating errors...');
    my %build_errors;
    my %guessed_errors;
    my $models_with_errors = 0;
    my $models_with_guessed_errors = 0;
    for my $build ( values %builds ) {
        my $key = 'Unknown';
        my $msg = 'Failure undetermined!';
        my $error = $self->_pick_optimal_error_log($build);
        if ( $error
                and
            ( ($error->file and $error->line) or ($error->inferred_file and $error->inferred_line) )
                and
            ($error->message or $error->inferred_message)
        ) {
            if ( $error->file and $error->line ) {
                $key = $error->file.' '.$error->line;
            } elsif ( $error->inferred_file and $error->inferred_line ) {
                $key = $error->inferred_file.' '.$error->inferred_line;
            } else {
                $key = 'unknown';
            }

            if ( $error->message ) {
                $msg = $error->message;
            } elsif ( $error->inferred_message ) {
                $msg = $error->inferred_message;
            } else {
                $msg = 'unknown';
            }

            $models_with_errors++;
        }
        elsif ( my $guessed_error = $self->_guess_build_error($build) ) {
            if ( not $guessed_errors{$guessed_error} ) {
                $guessed_errors{$guessed_error} = scalar(keys %guessed_errors) + 1;
            }
            $key = "Unknown, best guess #".$guessed_errors{$guessed_error};
            $msg = $guessed_error;
            $models_with_guessed_errors++;
        }
        $build_errors{$key} = "File:\n$key\nExample error:\n$msg\nModel\t\tBuild\t\tType/Failed Stage:\n" if not $build_errors{$key};
        my $type_name = $build->type_name;
        $type_name =~ s/\s+/\-/g;
        my %failed_events = map { $_->event_type => 1 } grep { $_->event_type ne 'genome model build' } $build->events('event_status in' => [qw/ Crashed Failed /]);
        my $failed_event = (keys(%failed_events))[0] || '';
        $failed_event =~ s/genome model build $type_name //;
        if ($failed_event eq '') {
            if ($build->status eq 'Unstartable') {
                $failed_event = 'Unstartable';
            }
        }
        $build_errors{$key} .= join("\t", $build->model_id, $build->id, $type_name.' '.$failed_event)."\n";
    }

    # Report
    my $models_in_tickets = map { @{$tickets{$_}} }keys %tickets;
    my $models_not_in_tickets = keys %builds;
    $self->status_message('Models: '.($models_in_tickets+ $models_not_in_tickets));
    $self->status_message('Models in tickets: '.$models_in_tickets);
    $self->status_message('Models not in tickets: '.$models_not_in_tickets);
    $self->status_message('Models with error log: '.$models_with_errors);
    $self->status_message('Models with guessed errors: '.$models_with_guessed_errors);
    $self->status_message('Models with unknown failures: '.($models_not_in_tickets - $models_with_errors - $models_with_guessed_errors));
    $self->status_message('Summarized errors: ');
    $self->status_message(join("\n", map { $build_errors{$_} } sort keys %build_errors));

    return 1;
}

sub get_events {
    my $self = shift;

    # Find cron models by failed build events
    my @events;
    if ($self->include_failed) {
        $self->status_message('Looking for failed model events...');
        @events = Genome::Model::Event->get(
            "event_status in" => ['Failed', 'Unstartable'],
            event_type => 'genome model build',
            user_name => 'apipe-builder',
            -hint => [qw/ build /],
        );
    }

    # Find cron models by unstartable build events
    if ($self->include_unstartable) {
        $self->status_message('Looking for unstartable model events...');
        my @unstartable_events = Genome::Model::Event->get(
            event_status => 'Unstartable',
            event_type => 'genome model build',
            user_name => 'apipe-builder',
            -hint => [qw/ build /],
        );
        @events = (@events, @unstartable_events);
    }

    if (scalar(@events)) {
        return @events;
    } else {
        $self->error_message('No failed or unstartable events found!');
        die $self->error_message();
    }
}

sub get_builds {
    my ($self, @events) = @_;

    $self->status_message("Looking up Models and Builds from events...");

    my %builds;
    for my $event (@events) {
        next if not $event->build_id;

        my $build = Genome::Model::Build->get(id => $event->build_id, -hint => [qw/ model events /]);
        my $model = $build->model;

        #If the latest build of the model succeeds, ignore those old
        #failing ones that will be cleaned by admin "cleanup-succeeded".
        if ($model->status) {
            next if $model->status eq 'Succeeded';

            unless ($self->include_pending) {
                next if $model->status eq 'Requested';
                next if $model->status eq 'Scheduled';
                next if $model->status eq 'Running';
            }
        }

        # only keep the most recently scheduled build
        next if $builds{ $model->id } and $builds{ $model->id }->date_scheduled gt $build->date_scheduled;
        $builds{ $model->id } = $build;
    }
    $self->status_message('Found '.keys(%builds).' models');
    return %builds;
}

sub get_tickets {
    my $self = shift;
    my %builds = %{(shift)};

    # Connect
    my $rt = _login_sso();

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
            query => "Queue = 'apipe-builds' AND ( Status = 'new' OR Status = 'open' )",

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
    $self->status_message('Found '.@ticket_ids.' tickets');

    # Go through tickets
    my %tickets;

    # re-set the login cookies that we saved away eariler
    $rt->_ua->cookie_jar($login_cookies);
    $self->status_message('Matching models and builds to tickets...');
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

        my $transactions = $ticket->transactions;
        my $transaction_iterator = $transactions->get_iterator;
        while ( my $transaction = &$transaction_iterator ) {
            my $content = $transaction->content;
            for my $model_id ( keys %builds ) {
                my $build_id = $builds{$model_id}->id;
                next if $content !~ /$model_id/ and $content !~ /$build_id/;
                push @{$tickets{$ticket_id.' '.$ticket->subject}}, $model_id;
            }
        }
    }

    $self->status_message('Found '.scalar(values(%tickets)).' models mentioned in the tickets');
    return %tickets;
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
    $mech->get( _server() );

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

sub _guess_build_error {
    my ($self, $build) = @_;

    if ($build->status eq 'Unstartable') {
        return $self->_guess_build_error_from_notes($build);
    }
    else {
        return $self->_guess_build_error_from_logs($build);
    }
}

sub _guess_build_error_from_notes {
    my ($self, $build) = @_;

    my $error = "Unstartable unknown error";
    my @notes = $build->notes;
    if (not @notes) {
        return $error;
    }
    my @unstartable_notes = grep {$_->header_text eq 'Unstartable'} @notes;
    if (not @unstartable_notes) {
        return $error;
    }
    my $note = $unstartable_notes[0];

    my $text = $note->body_text;
    if (not $text) {
        return $error;
    }

    my @lines = split(/\n/, $text);

    @lines = grep (/ERROR:\s+/, @lines);

    my $line_count = scalar @lines;
    if ($line_count == 0) {
        return $error;
    }

    my %errors;
    foreach my $line (@lines) {
        my ($err) = (split(/ERROR:\s+/, $line))[1];
        chomp $err;
        $errors{$err} = 1;
    }
    return join("\n", sort keys %errors);
}

sub _guess_build_error_from_logs {
    my ($self, $build) = @_;

    my $data_directory = $build->data_directory;
    my $log_directory = $data_directory.'/logs';
    return unless -d $log_directory;
    my %errors;
    find(
        sub{
            return unless $_ =~ /\.err$/;
            my @grep = (fgrep { /ERROR:\s+/ } $_ );
            return if $grep[0]->{count} == 0;
            for my $line ( values %{$grep[0]->{matches}} ) {
                my ($err) = (split(/ERROR:\s+/, $line))[1];
                chomp $err;
                next if $err eq "Can't convert workflow errors to build errors";
                next if $err eq 'relation "error_log_entry" does not exist';
                next if $err =~ /current transaction is aborted/;
                next if $err =~ /run_workflow_ls/;
                $errors{$err} = 1;
            }
        },
        $log_directory,
    );

    return join("\n", sort keys %errors);
}

sub _pick_optimal_error_log{
    my $self = shift;
    my $build = shift;
    my @errors = Genome::Model::Build::ErrorLogEntry->get(build_id => $build->id);
    my @optimal_errors = grep($_->file, @errors);
    unless (@optimal_errors){
        @optimal_errors = grep($_->inferred_file, @errors);
    }
    unless(@optimal_errors){
        return 0;
    }
    return shift @optimal_errors;
}
1;
