package Genome::Disk::Command::Allocation::UnarchiveBase;

use Genome;
use JSON 'to_json';

use strict;
use warnings;

class Genome::Disk::Command::Allocation::UnarchiveBase {
    is => 'Command::V2',
    is_abstract => 1,
    has => {
        lab => {
            is => 'Text',
            valid_values => [qw(Ding Informatics Maher Mardis-Wilson Mitreva Warren Weinstock)],
            doc => 'The lab which requests the data be unarchived.',
        },
        requestor => {
            is => 'Genome::Sys::User',
            is_optional => 1,
            doc => 'The user who is requesting the unarchive.',
        },
    },
};

sub execute {
    my $self = shift;

    unless ($self->requestor) {
        $self->requestor(Genome::Sys->current_user);
    }
    $self->status_message(sprintf("Unarchiving allocation(s) for the %s lab at the request of %s.",
            $self->lab, $self->requestor->name));

    return $self->_execute();
}

sub reason {
    my $self = shift;

    my %reason = (
        'lab'=>$self->lab,
        'requestor'=> $self->requestor->id
    );
    return to_json(\%reason);
}

1;
