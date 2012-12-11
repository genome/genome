package Genome::InstrumentData::Command::Imported::BlacklistSegments;

use strict;
use warnings;

use Genome;
use Set::Scalar;
use Carp "confess";


class Genome::InstrumentData::Command::Imported::BlacklistSegments {
    is => 'Command::V2',

    has_input => [
        imported_instrument_data => {
            is => 'Genome::InstrumentData::Imported',
            shell_args_position => 1,
            doc => "InstrumentData which has segments you'd like to blacklist."
        },
        segments => {
            is => 'Text',
            is_many => 1,
            shell_args_position => 2,
            doc => 'the name of segments to be blacklisted or removed from ' .
                   'the blacklist if they are already blacklisted.',
        },
    ],
    doc => "blacklist segments of an imported instrument-data",
};

sub help_detail {
    return <<EOS;
Add/Remove segments from the blacklist of segments.  Blacklisted
segments are not returned when you call get_segments() so are
not aligned in ReferenceAlignment builds.  All
segments can be found by passing allow_blacklisted_segments => 1
to get_segments().
EOS
}

sub execute {
    my $self = shift;

    my $instrument_data = $self->imported_instrument_data;
    my @segments = $self->segments;
    my $segments = Set::Scalar->new(@segments);

    my $available_segments = $instrument_data->get_read_groups_set();
    my $already_blacklisted = Set::Scalar->new(
            $instrument_data->blacklisted_segments);

    $self->_validate_segments($segments, $available_segments);

    my $to_be_added   = $segments - $already_blacklisted; # difference
    my $to_be_removed = $segments * $already_blacklisted; # intersection

    my @to_be_added = $to_be_added->members();
    if (@to_be_added) {
        $self->status_message(sprintf("Adding %s to the blacklist.",
                join(", ", map {"$_"} @to_be_added)));
        for my $segment (@to_be_added) {
            $instrument_data->add_blacklisted_segment($segment);
        }
    }

    my @to_be_removed = $to_be_removed->members();
    if (@to_be_removed) {
        $self->status_message(sprintf("Removing %s from the blacklist.",
                join(", ", map {"$_"} @to_be_removed)));
        for my $segment (@to_be_removed) {
            $instrument_data->remove_blacklisted_segment($segment);
        }
    }

    my @blacklisted_segments = $instrument_data->blacklisted_segments;
    if (@blacklisted_segments) {
        $self->status_message("The blacklist is now:\n" .
                join(", ", @blacklisted_segments));
    } else {
        $self->status_message("The blacklist is now empty.");
    }

    return 1;
}

sub _validate_segments {
    my ($self, $segments, $available_segments) = @_;

    my $invalid_segments = $segments - $available_segments;
    my @invalid_segments = $invalid_segments->members();
    if (@invalid_segments) {
        my $error = sprintf("Found %d invalid segment(s):\n    %s\n" .
                "not in the set:\n    %s", scalar(@invalid_segments),
                join(", ", map {"$_"} @invalid_segments),
                join(", ", map {"$_"} $available_segments->members()));
        confess $error;
    }
}


1;
