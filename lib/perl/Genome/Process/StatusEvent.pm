package Genome::Process::StatusEvent;

use strict;
use warnings FATAL => 'all';
use Genome;
use Data::Dump qw(pp);

my $VALID_STATUS_TRANSITIONS = {
    New => Set::Scalar->new('Scheduled', 'Running'),
    Scheduled => Set::Scalar->new('Running'),
    Running => Set::Scalar->new('Crashed', 'Succeeded'),
    Crashed => Set::Scalar->new(),
    Succeeded => Set::Scalar->new(),
};
my $VALID_STATUS_VALUES = [keys %$VALID_STATUS_TRANSITIONS];

class Genome::Process::StatusEvent {
    table_name => 'process.status_event',
    schema_name => 'GMSchema',
    data_source => 'Genome::DataSource::GMSchema',
    id_by => [
        id => {is => 'Text'},
    ],
    has => [
        process => {
            is => 'Genome::Process',
            id_by => 'process_id',
        },
        process_id => {
            is => 'Text',
        },
        old_status => {
            is => 'Text',
            is_optional => 1,
            valid_values => $Genome::Process::VALID_STATUS_VALUES,
            doc => 'The status of the process before the event',
        },
        new_status => {
            is => 'Text',
            valid_values => $Genome::Process::VALID_STATUS_VALUES,
            doc => 'The status of the process after the event',
        },
        timestamp => {
            is => 'Timestamp',
            doc => "When the status changed from 'old_status' to 'new_status'",
        },
    ],
};

sub valid_transitions {
    my $self = shift;

    if (defined($self->old_status)) {
        return $VALID_STATUS_TRANSITIONS->{$self->old_status};
    } else {
        return Set::Scalar->new('New');
    }
}

sub create {
    my $class = shift;

    my $self = $class->SUPER::create(@_);
    return unless $self;

    unless ($self->valid_transitions->contains($self->new_status)) {
        $self->error_message("Invalid transition from (%s) to (%s) only %s " .
            "is allowed.",
            $self->old_status || 'undef',
            $self->new_status,
            pp($self->valid_transitions->members));
        return;
    }
    $self->timestamp(UR::Context->current->now);

    return $self;
}

1;
