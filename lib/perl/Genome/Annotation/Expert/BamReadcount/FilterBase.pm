package Genome::Annotation::Expert::BamReadcount::FilterBase;

use strict;
use warnings;
use Genome;
use Genome::File::BamReadcount::Entry;

class Genome::Annotation::Expert::BamReadcount::FilterBase {
    is => ['Genome::Annotation::Filter::Base', 'Genome::Annotation::Filter::WithSampleName'],
    has => [
    ],
};

sub get_readcount_entry {
    my $self = shift;
    my $entry = shift;

    my $bam_readcount_string = $entry->sample_field($self->sample_index($entry->{header}), 'BRCT');
    return unless $bam_readcount_string;
    return Genome::File::BamReadcount::Entry->new(
        Genome::File::BamReadcount::Entry::decode($bam_readcount_string));
}

1;

