package Genome::VariantReporting::Suite::BamReadcount::VafCalculator;

use strict;
use warnings;
use Genome;
use List::Util 'sum';

sub calculate_vaf_for_all_alts {
    my ($entry, $readcount_entries) = @_;

    my $alt_alleles = $entry->{alternate_alleles};
    my $ref = $entry->{reference_allele};

    my %vafs;
    for my $allele (@$alt_alleles) {
        if (defined $readcount_entries->{$allele}) {
            $vafs{$allele} = calculate_vaf($readcount_entries->{$allele}, $allele, $entry->{reference_allele});
        }
        else {
            $vafs{$allele} = undef;
        }
    }
    return %vafs;
}

sub calculate_per_library_vaf_for_all_alts {
    my ($entry, $readcount_entry) = @_;
    my $alt_alleles = $entry->{alternate_alleles};
    my $ref = $entry->{reference_allele};

    my $vafs;
    for my $allele (@$alt_alleles) {
        my $vaf = calculate_per_library_vaf($readcount_entry, $allele, $entry->{reference_allele});
        $vafs->{$allele} = $vaf;
    }
    return $vafs;
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
    my $coverage = calculate_coverage_for_allele($bam_readcount_entry, $alt_allele, $ref);
    my $depth = $bam_readcount_entry->depth;
    if ($depth == 0) {
        return 0;
    } else {
        return $coverage / $depth * 100;
    }
}

sub calculate_per_library_vaf {
    my ($bam_readcount_entry, $alt_allele, $ref) = @_;
    my $coverage = calculate_per_library_coverage_for_allele(@_);
    my $vaf;
    while ( my ($library, $readcount) = each %$coverage ) {
        my $library_depth = depth_for_library_name($bam_readcount_entry, $library);
        if ($library_depth == 0) {
            $vaf->{$library} = 0;
        }
        else {
            $vaf->{$library} = $readcount / $library_depth * 100;
        }
    }
    return $vaf;
}

sub depth_for_library_name {
    my ($bam_readcount_entry, $library_name) = @_;

    for my $library ($bam_readcount_entry->libraries) {
        if ($library->name eq $library_name) {
            return $library->depth;
        }
    }
    die sprintf("No library information for library named (%s)", $library_name);
}

sub calculate_coverage_for_allele {
    my $counts = calculate_per_library_coverage_for_allele(@_);
    return sum(values %$counts);
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

sub translated_allele {
    my ($allele, $ref) = @_;

    my $translated_allele;
    if (is_insertion($ref, $allele)) {
        $translated_allele = "+".translate_pure_indels($allele, $ref);
    }
    elsif (is_deletion($ref, $allele)) {
        $translated_allele = "-".translate_pure_indels($ref, $allele);
    }
    else {
        $translated_allele = $allele;
    }
    return $translated_allele;
}

#This code attempts to calculate which bases of the longer sequence were
#removed to arrive at the shorter sequence.
#For a deletion the $long sequence is the reference and the $short sequence is
#the alternate allele. We return the bases that were removed from the
#reference to arrive at the alternate allele.
#For an insertion the $short sequence is the reference and the $long sequence is
#the alternate allele. We return the bases that were added to the reference
#to arrive at the alternate allele
sub translate_pure_indels {
    my ($long, $short) = @_;
    my @long = split('', $long);
    my @short = split('', $short);

    #This will match up the $short sequence to the $long sequence, starting
    #from the end. It will then calculate the first position (inclusive)
    #that differs. This is the end position of the indel.
    #We remove all matched bases from the $short sequence. This leaves only
    #bases that still need to be matched the $long sequence, starting from the
    #beginning.
    my $end = $#long;
    for my $i (reverse 0..$#long) {
        if ($long[$i] eq $short[$#short]) {
            $end = $i - 1;
            pop @short;
            last if (scalar(@short) == 0);
        }
        else {
            last;
        }
    }

    #This will match up the $short sequence to the $long sequence, starting
    #from the beginning. It will then calculate the first position (inclusive)
    #that differs. This is the start position of the indel.
    my $start = 0;
    for my $i (0..$#short) {
        if ($long[$i] eq $short[$i]) {
            $start = $i + 1;
        }
        else {
            last;
        }
    }
    return substr($long, $start, $end - $start + 1);
}

1;
