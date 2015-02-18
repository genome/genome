package Genome::Notify::Notifier::Ignore;

use strict;
use warnings;

use Genome;

class Genome::Notify::Notifier::Ignore {
    is => 'Genome::Notify::Notifier',
    has_classwide_constant => [
        name => {
            value => 'ignore',
        },
   ],
   doc => 'The notice will be immediately acknowledged.',
};

sub notify {
    my $class = shift;
    my $notice = shift;

    $notice->acknowledged(1);
}

1;
