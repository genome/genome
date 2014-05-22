package Genome::Annotation::Expert::BamReadcount::VafCalculator;

use strict;
use warnings;
use Genome;

sub calculate_vaf {
    my ($bam_readcount_entry, $alt_allele) = @_;

    my $var_count = 0;
    for my $lib ($bam_readcount_entry->libraries) {
        $var_count += $lib->metrics_for($alt_allele)->count;
    }

    return $var_count / $bam_readcount_entry->depth * 100;
}

1;

