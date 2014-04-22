package Genome::Annotation::BamReadcount::MinCoverageFilter;

use strict;
use warnings FATAL => 'all';
use Genome;
use List::Util qw/first/;

class Genome::Annotation::BamReadcount::MinCoverageFilter {
    is => 'Genome::Annotation::FilterBase',
    has => {
    },
};

sub name {
    return 'min-coverage';
}

sub requires_experts {
    return ();
}

sub process_entry {
    my ($self, $entry) = @_;

    return 0;
}

1;
