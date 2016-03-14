package Genome::Model::Command::Admin::FailedModelTickets;

use strict;
use warnings;

use Genome;

BEGIN {
        $ENV{UR_DBI_NO_COMMIT} = 1;
}

class Genome::Model::Command::Admin::FailedModelTickets {
    is => 'Genome::Model::Command::Admin::FailedModelTicketBase',
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
This command collects cron models with failed or unstartable builds and scours tickets for their IDs. If they are not found in the existing tickets, the models are summarized and grouped by the `genome model build determine-error` output.
HELP
}

sub execute {
    my $self = shift;

    my %builds = $self->get_builds();

    $self->remove_builds_in_tickets(\%builds);

    my $build_errors = $self->get_build_errors(%builds);

    # Report
    $self->status_message("\n\n");
    $self->status_message(join("\n\n", map { $build_errors->{$_} } sort keys %$build_errors));

    return 1;
}

sub get_builds {
    my $self = shift;

    my $production_role = Genome::Sys::User::Role->get(name => 'production');
    my @prod_users = $production_role->users;
    my @usernames = map $_->username, @prod_users;

    # Find cron models by failed build events
    my @builds;
    if ($self->include_failed) {
        $self->status_message('Looking for failed builds...');
        @builds = Genome::Model::Build->get(
            status => 'Failed',
            'model.run_as' => \@usernames,
        );
    }

    # Find cron models by unstartable build events
    if ($self->include_unstartable) {
        $self->status_message('Looking for unstartable builds...');
        my @unstartable_builds = Genome::Model::Build->get(
            status => 'Unstartable',
            'model.run_as' => \@usernames,
        );
        @builds = (@builds, @unstartable_builds);
    }

    unless (scalar(@builds)) {
        $self->error_message('No failed or unstartable builds found!');
        die $self->error_message();
    }

    my @models = Genome::Model->get('id in' => [map {$_->model_id} @builds]);
    $self->status_message(sprintf("Found %d builds in %d models",
            scalar(@builds), scalar(@models)));

    # cache the model_status calculations (they're calculated and slowly)...
    my %model_status;
    for my $model (@models) {
        $model_status{$model->id} = $model->status;
    }

    $self->status_message("Filtering down to latest build for each model...");
    my %builds;
    for my $build (@builds) {
        my $model = $build->model;

        #If the latest build of the model succeeds, ignore those old
        #failing ones that will be cleaned by admin "cleanup-succeeded".
        my $model_status = $model_status{$model->id};
        if ($model_status) {
            next if $model_status eq 'Succeeded';

            unless ($self->include_pending) {
                next if $model_status eq 'Requested';
                next if $model_status eq 'Scheduled';
                next if $model_status eq 'Running';
            }
        }

        # only keep the most recently scheduled build
        next if $builds{ $model->id } and $builds{ $model->id }->created_at gt $build->created_at;
        $builds{ $model->id } = $build;
    }
    $self->status_message('Found '.keys(%builds).' models');
    return %builds;
}


sub get_build_errors {
    my ($self, %builds) = @_;

    $self->status_message('Categorizing builds...');

    my %build_errors;

    for my $build ( values %builds ) {
        my $cmd = Genome::Model::Build::Command::DetermineError->execute(
            build => $build,
            display_results => 0,
        );

        my $key = "Unknown Error";
        my $header = join("\t", qw(Model Build Build-Class Date));
        my $line = join("\t", $build->model_id, $build->id, $build->class, $cmd->error_date);
        if ($cmd->error_type eq 'Unstartable') {
            $key = $cmd->get_unstartable_key;
        } elsif ($cmd->error_type eq 'Failed') {
            $key = $cmd->get_failed_key;

            # unstartable and unknown errors don't generally have a host.
            $header = join("\t", qw(Model Build Build-Class Host Date));
            $line = join("\t", $build->model_id, $build->id, $build->class, $cmd->error_host, $cmd->error_date);
        }

        unless ($build_errors{$key}) {
            $build_errors{$key} = sprintf("%s: %s\n%s:\n%s\n\n%s\n",
                $self->_color("Key:", "bold"), $key,
                $self->_color("Example error:", "bold"), $cmd->error_text,
                $header,
            );
        }
        $build_errors{$key} .= $line . "\n";
    }

    return \%build_errors;
}

sub remove_builds_in_tickets {
    my ($self, $builds) = @_;

    # Connect
    my $rt = $self->_login_sso();

    # Go through tickets
    my @ticket_ids = $self->_find_open_tickets($rt);
    my %tickets;

    $self->status_message('Matching models and builds to tickets...');
    for my $ticket_id ( @ticket_ids ) {
        my $ticket = $self->_ticket_for_id($rt, $ticket_id);
        next unless $ticket;

        my $transactions = $ticket->transactions;
        my $transaction_iterator = $transactions->get_iterator;
        while ( my $transaction = &$transaction_iterator ) {
            my $content = $transaction->content;
            next unless $content;
            for my $model_id ( keys %$builds ) {
                my $build_id = $builds->{$model_id}->id;
                next if $content !~ /$model_id/ and $content !~ /$build_id/;
                push @{$tickets{$ticket_id.' '.$ticket->subject}}, $model_id;
            }
        }
    }

    for my $model_id_list (values %tickets) {
        for my $model_id (@$model_id_list) {
            delete $builds->{$model_id};
        }
    }

    my $models_in_tickets = map { @{$tickets{$_}} } keys %tickets;
    my $models_not_in_tickets = keys %$builds;
    $self->status_message($self->_color('Tickets mentioning models/builds: ', 'bold').scalar(values(%tickets)));
    $self->status_message($self->_color('Models: ', 'bold').($models_in_tickets + $models_not_in_tickets));
    $self->status_message($self->_color('Models in tickets: ', 'bold').$models_in_tickets);
    $self->status_message($self->_color('Models not in tickets: ', 'bold').$models_not_in_tickets);

    return %tickets;
}


1;
