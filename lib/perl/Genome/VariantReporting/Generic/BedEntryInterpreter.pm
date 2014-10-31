package Genome::VariantReporting::Generic::BedEntryInterpreter;

use strict;
use warnings;
use Genome;
use Genome::Utility::Vcf qw(convert_indel_gt_to_bed);

class Genome::VariantReporting::Generic::BedEntryInterpreter {
    is => ['Genome::VariantReporting::Framework::Component::Interpreter'],
};

sub name {
    return 'bed-entry';
}

sub requires_annotations {
    ();
}

sub field_descriptions {
    return (
        chromosome_name => 'Chromosome: an identifier from the reference genome or an angle-bracketed ID String ("<ID>") pointing to a contig in the assembly file',
        start => 'The start position of the variant. Zero-based',
        stop => 'The end position of the variant. Zero-based, exclusive',
    );
}

sub _interpret_entry {
    my ($self, $entry, $passed_alt_alleles) = @_;

    my %return_values;
    for my $alt_allele (@$passed_alt_alleles) {
        my ($chromosome, $start, $stop) = $self->convert_to_bed($entry, $alt_allele);
        $return_values{$alt_allele} =  {
            chromosome_name => $chromosome,
            start           => $start,
            stop            => $stop,
        };
    }
    return %return_values;
}

sub convert_to_bed {
    my $self = shift;
    my $entry = shift;
    my $variant_allele = shift;

    if (length($variant_allele) == length($entry->{reference_allele})) {
        return ($entry->{chrom}, $entry->{position}-1, $entry->{position});
    }
    else {
        my ($indel_gts, $indel_shifts) = convert_indel_gt_to_bed($entry->{reference_allele}, $variant_allele);
        unless (scalar(@$indel_gts) == 1 && scalar(@$indel_shifts) == 1) {
            die $self->error_message("Expected one genotype and one shift. Got (%s) and (%s)", scalar(@$indel_gts), scalar(@$indel_shifts));
        }
        my $indel = $indel_gts->[0];
        my ($ref, $var) = @$indel;
        my $shift = $indel_shifts->[0];
        my $pos = $entry->{position} + $shift;
        if($ref eq '*') {
            #insertion
            return ($entry->{chrom}, $pos-1, $pos-1);
        }
        if($var eq '*') {
            #deletion
            return ($entry->{chrom}, $pos-1, $pos-1 + length($ref));
        }
    }
}

1;
