package Genome::Process::Command::ListRecent;

use strict;
use warnings FATAL => 'all';
use Genome;
use Time::Piece qw();
use Time::Seconds;

class Genome::Process::Command::ListRecent {
    is => 'Genome::Object::Command::List',
    has => [
        subject_class_name  => {
            is_constant => 1,
            value => 'Genome::Process'
        },
        user => {
            is => 'Genome::Sys::User',
            is_optional => 1,
            shell_args_position => 1,
            doc => 'The user to list processes for',
        },
        days => {
            is => 'Number',
            is_optional => 1,
            default_value => 3,
            shell_args_position => 2,
            doc => 'How many days in the past to look for ended processes',
        },
        show => { default_value => 'id,status,created_by,metadata_directory' },
    ],
    has_optional_transient => [
        filter => {
            is => 'Text',
        },
        show => {
            is => 'Text',
        },
    ],
    doc => 'list genome processes',
};

my $date_pattern ='%Y-%M-%d %H:%M:%S%Y';

sub execute{
    my $self = shift;

    $self->status_message("\nProcesses that finished in the last %0.2f days:",
        $self->days,
    );
    $self->show('id,subclass_name,status,ended_at,metadata_directory');
    $self->order_by('ended_at');

    my $user = $self->get_user();
    $self->filter(join(',',
        sprintf('created_by=%s', $user->username),
        "status:Crashed/Succeeded",
        sprintf("ended_at>%s", $self->get_since_date),
    ));
    $self->SUPER::_execute_body(@_);


    $self->status_message("\nActive processes:");
    $self->show('id,subclass_name,status,started_at,metadata_directory');
    $self->order_by('started_at');
    $self->filter(join(',',
        sprintf('created_by=%s', $user->username),
        "status:New/Scheduled/Running",
    ));
    $self->SUPER::_execute_body(@_);

    return 1;
}

sub get_user {
    my $self = shift;
    if (defined($self->user)) {
        return $self->user;
    } else {
        return Genome::Sys->current_user;
    }
}

sub get_since_date {
    my $self = shift;

    return Time::Piece->localtime() - ONE_DAY*$self->days;
}


1;
