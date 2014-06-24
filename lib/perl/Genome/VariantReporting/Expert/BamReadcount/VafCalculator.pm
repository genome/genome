package Genome::VariantReporting::Expert::BamReadcount::VafCalculator;

use strict;
use warnings;
use Genome;

sub calculate_vaf_for_all_alts {
    my ($entry, $readcount_entry) = @_;
    my $alt_alleles = $entry->{alternate_alleles};
    my $ref = $entry->{reference_allele};
    my @alleles;
    for my $allele (@$alt_alleles) {
        if (is_insertion($ref, $allele)) {
            push @alleles, "+".substr($allele, length($ref));
        }
        else {
            push @alleles, $allele;
        }
    }
    my %alleles = calculate_vaf_for_multiple_alleles($readcount_entry, \@alleles);
    my %translated_alleles;
    for my $allele (keys %alleles) {
        my $new_allele = $allele;
        $new_allele =~ s/^[+-]/$ref/;
        $translated_alleles{$new_allele} = $alleles{$allele};
    }
    return %translated_alleles;
}

sub is_insertion {
    my ($ref, $allele) = @_;
    if (length($ref) < length($allele)) {
        return 1;
    }
    return 0;
}

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

    return calculate_coverage_for_allele($bam_readcount_entry, $alt_allele) / $bam_readcount_entry->depth * 100;
}

sub calculate_coverage_for_allele {
    my ($bam_readcount_entry, $allele) = @_;

    my $count = 0;
    for my $lib ($bam_readcount_entry->libraries) {
        my $metrics;
        eval {
            $metrics = $lib->metrics_for($allele);
        };
        if (defined $metrics) {
            $count += $metrics->count;
        }
    }

    return $count;
}

1;

