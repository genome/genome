package Genome::File::Vcf::BamReadcountUtilities;

use strict;
use warnings FATAL => 'all';
use Genome;
use Genome::Utility::Vcf qw(convert_indel_gt_to_bed);

sub vcf_entry_to_allele_offsets {
    my ($entry) = @_;
    my @offsets;
    for my $allele (@{$entry->{alternate_alleles}}) {
        if (length($allele) == length($entry->{reference_allele})) {
           push @offsets, 0;
        }
        else {
            my (undef, $shifts) = convert_indel_gt_to_bed($entry->{reference_allele}, $allele);
            if (defined $shifts->[0]) { 
                if ($entry->is_deletion($allele)) {
                    push @offsets, $shifts->[0];
                }
                else {
                    push @offsets, $shifts->[0] - 1;
                }
            }
            else {
                push @offsets, 0;
            }
        }
    }
    return @offsets;
}


sub vcf_entry_to_readcount_positions {
    my ($entry) = @_;
    my $pos = $entry->{position};
    my %positions = map { $pos + $_ => 1 } vcf_entry_to_allele_offsets($entry);
    return sort { $a <=> $b } keys %positions;
}

sub entries_match {
    my ($readcount_entry, $vcf_entry) = @_;

    if (defined $readcount_entry) {
        my $rc_chrom  = $readcount_entry->chromosome;
        my $rc_pos    = $readcount_entry->position;
        my $vcf_chrom = $vcf_entry->{chrom};
        if ($rc_chrom eq $vcf_chrom) {
            my @positions = vcf_entry_to_readcount_positions($vcf_entry);
            for my $vcf_pos (@positions) {
                if ($rc_pos == $vcf_pos) {
                    return 1;
                }
            }
        }
    }
    return 0;
}

1;
