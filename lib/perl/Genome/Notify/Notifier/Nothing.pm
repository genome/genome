package Genome::Notify::Notifier::Nothing;

use strict;
use warnings;

use Genome;

class Genome::Notify::Notifier::Nothing {
    is => 'Genome::Notify::Notifier',
    has_classwide_constant => [
        name => {
            value => 'do nothing',
        },
   ],
   doc => 'The notice will not be actively communicated, but will remain unacknowledged.',
};

sub notify {
    my $class = shift;
    my $notice = shift;

    return 1;
}

1;
