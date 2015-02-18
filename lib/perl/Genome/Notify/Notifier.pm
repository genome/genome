package Genome::Notify::Notifier;

use strict;
use warnings;

use Genome;

class Genome::Notify::Notifier {
    is => 'UR::Object',
    is_abstract => 1,
    has_classwide => [
        name => {
            is => 'Text',
            doc => 'A user-friendly name for this notifier',
        },
    ],
};

sub notify {
    die 'subclass must implement notify';
}

1;
