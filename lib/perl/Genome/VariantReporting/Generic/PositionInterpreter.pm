package Genome::VariantReporting::Generic::PositionInterpreter;

use strict;
use warnings;
use Genome;
use Genome::Utility::Vcf qw(convert_indel_gt_to_bed);

class Genome::VariantReporting::Generic::PositionInterpreter {
    is => 'Genome::VariantReporting::Framework::Component::Interpreter',
    doc => 'Output the chromosome name, the start position, the stop position, the reference allele, and the variant allele',
};

sub name {
    return 'position';
}

sub requires_annotations {
    ();
}

sub field_descriptions {
    return (
        chromosome_name => 'Chromosome: an identifier from the reference genome or an angle-bracketed ID String ("<ID>") pointing to a contig in the assembly file',
        start => 'The start position of the variant. One-based',
        stop => 'The end position of the variant. One-based, inclusive',
        reference => 'The reference allele at this position',
        variant => 'The variant called at this position',
    );
}

sub _interpret_entry {
    my $self = shift;
    my $entry = shift;
    my $passed_alt_alleles = shift;

    my %return_values;
    for my $variant_allele (@$passed_alt_alleles) {
        my ($start, $stop, $reference, $variant);
        if (length($variant_allele) != length($entry->{reference_allele})) {
            ($start, $stop, $reference, $variant) = $self->_interpret_indel($entry, $variant_allele);
        }
        else {
            $start = $entry->{position};
            $stop = $entry->{position} + length($variant_allele) - 1;
            $reference = $entry->{reference_allele};
            $variant = $variant_allele;
        }
        $return_values{$variant_allele} = {
            chromosome_name => $entry->{chrom},
            start           => $start,
            stop            => $stop,
            reference       => $reference,
            variant         => $variant,
        };
    }

    return %return_values;
}

sub _interpret_indel {
    my ($self, $entry, $variant_allele) = @_;
    my ($indel_gts, $indel_shifts) = convert_indel_gt_to_bed($entry->{reference_allele}, $variant_allele);
    if (@$indel_gts != 1) {
        die $self->error_message("Vcf indel did not translate to exactly one bed genotype");
    }
    my $indel = $indel_gts->[0];
    my ($ref, $var) = @$indel;
    my $shift = $indel_shifts->[0];

    my $start = $entry->{position} + $shift;
    my ($stop, $reference, $variant);
    if ($ref eq '*') {
        $ref = '-';
        $start--;
        $stop = $start+1;
    }
    elsif ($var eq '*') {
        $var = '-';
        $stop = $start + length($ref) - 1;
    }
    else {
        die $self->error_Message("Vcf indel translated with unexpected ref and var");
    }
    return ($start, $stop, $ref, $var);
}

1;
