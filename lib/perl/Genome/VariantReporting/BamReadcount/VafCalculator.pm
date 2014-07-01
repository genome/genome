package Genome::VariantReporting::BamReadcount::VafCalculator;

use strict;
use warnings;
use Genome;

sub calculate_vaf_for_all_alts {
    my ($entry, $readcount_entry) = @_;
    my $alt_alleles = $entry->{alternate_alleles};
    my $ref = $entry->{reference_allele};
    my @alleles;
    my %return_values;
    for my $allele (@$alt_alleles) {
        my $vaf = calculate_vaf($readcount_entry, $allele, $ref);
        $return_values{$allele} = $vaf;
    }
    return %return_values;
}

sub is_insertion {
    my ($ref, $allele) = @_;
    if (length($ref) < length($allele)) {
        return 1;
    }
    return 0;
}

sub is_deletion {
    my ($ref, $allele) = @_;
    if (length($ref) > length($allele)) {
        return 1;
    }
    return 0;
}

sub calculate_vaf {
    my ($bam_readcount_entry, $alt_allele, $ref) = @_;

    return calculate_coverage_for_allele($bam_readcount_entry, $alt_allele, $ref) / $bam_readcount_entry->depth * 100;
}

sub translated_allele {
    my ($allele, $ref) = @_;
    my $translated_allele;
    if (is_insertion($ref, $allele)) {
        $translated_allele = "+".substr($allele, length($ref));
    }
    elsif (is_deletion($ref, $allele)) {
        $translated_allele = "-".substr($ref, length($allele));
    }
    else {
        $translated_allele = $allele;
    }
    return $translated_allele;
}

sub calculate_coverage_for_allele {
    my ($bam_readcount_entry, $allele, $ref) = @_;

    my $count = 0;
    for my $lib ($bam_readcount_entry->libraries) {
        $count+=calculate_coverage_for_allele_and_library($allele, $ref, $lib);
    }

    return $count;
}

sub calculate_per_library_coverage_for_allele {
    my ($bam_readcount_entry, $allele, $ref) = @_;

    my $counts;
    for my $lib ($bam_readcount_entry->libraries) {
        $counts->{$lib->name} = calculate_coverage_for_allele_and_library($allele, $ref, $lib);
    }

    return $counts;
}

sub calculate_coverage_for_allele_and_library {
    my ($allele, $ref, $library) = @_;

    my $metrics;
    eval {
        $metrics = $library->metrics_for(translated_allele($allele, $ref));
    };
    if (defined $metrics) {
        return $metrics->count;
    } else {
        return 0;
    }
}

1;

