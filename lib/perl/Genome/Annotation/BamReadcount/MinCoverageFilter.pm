package Genome::Annotation::BamReadcount::MinCoverageFilter;

use strict;
use warnings FATAL => 'all';
use Genome;
use List::Util qw/first/;

class Genome::Annotation::BamReadcount::MinCoverageFilter {
    is => 'Genome::Annotation::FilterBase',
    has => {
        min_coverage => {
            is => 'Number',
        },
        sample_index => {
            is => 'Integer',
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

sub get_readcount_entry {
    my $self = shift;
    my $entry = shift;

    my $bam_readcount_string = $entry->sample_field($self->sample_index, 'BRCT');
    return unless $bam_readcount_string;
    return Genome::File::BamReadcount::Entry->new(
        Genome::File::BamReadcount::Entry::decode($bam_readcount_string));
}

1;
