package Genome::Model::Build::Command::AbandonFailed;

use strict;
use warnings;
use Genome;
use Carp;
use Benchmark;

class Genome::Model::Build::Command::AbandonFailed {
    is => 'Command',
    has => {
        user => {
            is => 'Text',
            is_optional => 1,
        },
    },
    doc => 'Abandons all failed builds by the user that are not the most recent buildfor a model',
};

sub help_brief {
    return 'Remove all abandoned builds owned by the user';
}

sub help_detail {
    return <<EOS
This command will grab all builds owned by the user and remove them all
EOS
}

sub execute {
    my $self = shift;

    my $alpha = Benchmark->new();
    my $user = $self->user;
    $user ||= Genome::Sys->username;
    confess "Could not get user name from environment" unless defined $user;

    my @builds = map {$_->build} Genome::Model::Event->get(
        event_type => 'genome model build',
        event_status => "Failed",
        user_name => $user,
        -hints => ["build"]
    );
    my $beta = Benchmark->new();
    #print  "got builds ".timestr(timediff($beta, $alpha))."\n";
    unless (@builds) {
        $self->status_message("No builds owned by $user are failed and eligible, so no removal necessary!");
        return 1;
    }
    my @eligible_to_abandon;
    for my $build (@builds){
        push @eligible_to_abandon, $build unless $build->id == $build->model->latest_build->id;
    }
    my $gamma = Benchmark->new();
    #print "pruned ineligible builds ".timestr(timediff($gamma,$beta))."\n";

    $self->status_message("User $user has ". scalar @builds ." failed builds.  Of these ". scalar @eligible_to_abandon ." are eligible to abandon.");
    $self->status_message(join("\n", map { $_->__display_name__ } @eligible_to_abandon));
    # TODO: use Genome::Command::Base->ask_user_question
    $self->status_message("Are you sure you want to abandon these builds? (y/n) [n]");
    my $answer = <STDIN>;
    return unless $answer =~ /^[yY]/;


    my $abandon_failed = Genome::Model::Build::Command::Abandon->create(
        builds => \@eligible_to_abandon,
    );

    return $abandon_failed->execute();
}

1;

