package Genome::Annotation::BamReadcount::MinCoverageFilter;

use strict;
use warnings FATAL => 'all';
use Genome;
use List::Util qw/first/;

class Genome::Annotation::BamReadcount::MinCoverageFilter {
    is => 'Genome::Annotation::BamReadcount::FilterBase',
    has => {
        min_coverage => {
            is => 'Number',
        },
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

    my $readcount_entry = $self->get_readcount_entry($entry);
    return 0 unless $readcount_entry;

    if ($readcount_entry->depth >= $self->min_coverage) {
        return 1;
    }
    return 0;
}

1;
