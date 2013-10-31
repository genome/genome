package Genome::Disk::Detail::StrictObject;

use strict;
use warnings;

use UR;
use Carp qw(confess);


class Genome::Disk::Detail::StrictObject {
    is_abstract => 1,
};


sub create {
    my $class = shift;

    my $self = $class->SUPER::create(@_);

    if ($self->__errors__) {
        my @messages;
        for my $error_tag ($self->__errors__) {
            push @messages, $error_tag->__display_name__ . "\n";
        }
        confess sprintf('Could not create object (%s):\n%s',
            $class, join('', @messages));
    }

    $self->sanitize;
    $self->validate;

    return $self;
}


sub sanitize {};
sub validate {};


1;
