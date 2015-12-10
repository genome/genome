package Genome::Disk::Command::Allocation::DelayArchiving;

use strict;
use warnings;
use Genome;
use DateTime;

class Genome::Disk::Command::Allocation::DelayArchiving {
    is => 'Command::V2',
    has => [
        allocations => {
            is => 'Genome::Disk::Allocation',
            is_many => 1,
            doc => 'allocations that are not to be archived, resolved via Command::V2',
        },
        _duration_property_for_commands(),
    ],
};

sub help_detail { 
    return 'Archiving of the given allocations and paths is delayed by 3 months (or the value you specified). Any paths that are given are resolved ' .
        'to allocations if possible, otherwise a warning is emitted.';
}
sub help_brief { return 'given allocations are marked so they cannot be archived' };
sub help_synopsis { return help_brief() . "\n" };

sub execute {
    my $self = shift;
    for my $allocation ($self->allocations) {
        $allocation->archive_after_time($self->_resolve_date_from_months($self->duration));
    }
    return 1;
}

sub _resolve_date_from_months {
    my ($class, $months) = @_;
    return DateTime->now(time_zone => 'local')
        ->add(months => $months)
        ->strftime('%F %T');
}

sub _duration_property_for_commands {
    (
        duration => {
            is => 'Number',
            default_value => 3,
            doc => 'The number of months to delay archiving',
            valid_values => [1..12],
        }
    )
}

1;

