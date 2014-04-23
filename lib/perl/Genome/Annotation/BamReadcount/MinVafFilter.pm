package Genome::Annotation::BamReadcount::MinVafFilter;

use strict;
use warnings;
use Genome;

class Genome::Annotation::BamReadcount::MinVafFilter {
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

    my %return_values;
    for my $alt_allele (@{$entry->{alternate_alleles}}) {
        $return_values{$alt_allele} = 0;
    }

    my @sample_alt_alleles = sort $entry->alt_bases_for_sample($self->sample_index);
    for my $sample_alt_allele (@sample_alt_alleles) {
        my $vaf = $self->calculate_vaf($entry, $sample_alt_allele);
        if ($vaf >= $self->min_vaf) {
            $return_values{$sample_alt_allele} = 1;
        }
    }

    return %return_values;
}

sub calculate_vaf {
    my $self = shift;
    my ($entry, $alt_allele) = @_;
    my $bam_readcount_entry = $self->get_readcount_entry($entry);

    my $var_count = 0;
    for my $lib ($bam_readcount_entry->libraries) {
        $var_count += $lib->metrics_for($alt_allele)->count;
    }

    return $var_count / $bam_readcount_entry->depth * 100;
}

1;

