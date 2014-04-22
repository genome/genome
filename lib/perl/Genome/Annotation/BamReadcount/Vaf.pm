package Genome::Annotation::BamReadcount::Vaf;

use strict;
use warnings;
use Genome;

class Genome::Annotation::BamReadcount::Vaf {
    is => 'Genome::Annotation::BamReadcount::FilterBase',
    has => [
        min_vaf => {
            is => 'Number',
        },
    ],
};

sub name {
    return 'min-vaf';
}

sub requires_experts {
    return ('bam-readcount');
}

sub process_entry {
    my $self = shift;
    my $entry = shift;

    my @alt_alleles = $entry->alt_bases_for_sample($self->sample_index);
    my $vaf = $self->calculate_vaf($entry, @alt_alleles);
    return $vaf >= $self->min_vaf;
}

sub calculate_vaf {
    my $self = shift;
    my ($entry, @alt_alleles) = @_;
    my $bam_readcount_entry = $self->get_readcount_entry($entry);

    my $var_count = 0;
    for my $lib ($bam_readcount_entry->libraries) {
        for my $allele (@alt_alleles) {
            $var_count += $lib->metrics_for($allele)->count;
        }
    }

    return $var_count / $bam_readcount_entry->depth * 100;
}

1;

