package Genome::ConfigValidatorBase;

use strict;
use warnings;

use Mouse::Role;

requires qw(
    check
    message
);

=item $validator->check($value)

Returns true iff the value passes the validator's constraint.

=item $validator->message()

Returns the message to display to a user if the C<check> fails.

=item $validator->validate($value)

Returns the error message for the value; returns undef if the value passes the validator's constraint.

=cut

sub validate {
    my ($self, $value) = @_;

    unless ($self->check($value)) {
        return $self->message;
    }

    return;
}

1;
