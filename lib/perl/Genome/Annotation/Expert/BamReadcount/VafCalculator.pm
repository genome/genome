package Genome::Annotation::Expert::BamReadcount::VafCalculator;

use strict;
use warnings;
use Genome;

sub calculate_vaf_for_multiple_alleles {
    my $bam_readcount_entry = shift;
    my $alt_alleles = shift;
    my %return_values;
    for my $sample_alt_allele (@$alt_alleles) {
        $return_values{$sample_alt_allele} = calculate_vaf(
            $bam_readcount_entry, $sample_alt_allele);
    }
    return %return_values;
}

sub calculate_vaf {
    my ($bam_readcount_entry, $alt_allele) = @_;

    my $var_count = 0;
    for my $lib ($bam_readcount_entry->libraries) {
        $var_count += $lib->metrics_for($alt_allele)->count;
    }

    return $var_count / $bam_readcount_entry->depth * 100;
}

1;

