package Genome::Model::Tools::Bed::Convert::CoordinateConverter;

use strict;
use warnings;
use Genome;
use Exporter 'import';

our @EXPORT_OK = qw(convert_indel_gt_to_bed);

sub convert_indel_gt_to_bed {
    my ($reference_allele, @alleles) = @_;

    #the challenge here is that VCF can represent multi-base calls in both the reference and variant and we clearly do NOT support that in our tgi-bed representation

    my %alleles = map { $_ => 1 } @alleles; #our bed doesn't care about indel genotypes.

    my @bed_alleles;
    my @position_shifts;
    for my $allele (keys %alleles) {
        next if $allele eq $reference_allele;
        my ($ref, $var, $right_shift) = _simplify_indel_allele($reference_allele, $allele);
        unless($var eq q{} || $ref eq q{}) {
            warn "Complex indels cannot be converted to TGI bed. This indel ($ref, $var) will be skipped.\n";
            next;
        }
        ($ref, $var) = map { $_ ne q{} ? $_ : '*' } ($ref, $var);

        push @bed_alleles, [$ref,$var];
        push @position_shifts, $right_shift;
    }
    return (\@bed_alleles, \@position_shifts);
}

sub _simplify_indel_allele {
    my ($ref, $var) = @_;
    #these could be padded e.g. G, GT for a T insertion or GCC G for a 2bp deletion
    #they could also be complex e.g. GCCCGT, GCGT for a 2bp deletion
    #they could also represent an insertion deletion event e.g. GCCCGT GCGGGGT; these cannot be represented in genome bed. Throw an error or warn.
    #
    #I think the algorithm should be trim end (no updating of coords required)
    #trim beginning and return number of bases trimmed

    my @ref_array = map { uc } split //, $ref;
    my @var_array = map { uc } split //, $var;

    while(@ref_array and @var_array and $ref_array[-1] eq $var_array[-1]) {
        pop @ref_array;
        pop @var_array;
    }

    my $right_shift = 0;
    while(@ref_array and @var_array and $ref_array[0] eq $var_array[0]) {
        shift @ref_array;
        shift @var_array;
        $right_shift++;
    }

    return (join("",@ref_array), join("",@var_array), $right_shift);
}

1;

