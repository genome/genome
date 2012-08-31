package Genome::Model::Tools::Vcf::Convert::Snv::VarscanSomatic;

use strict;
use warnings;
use Genome;
use Genome::Info::IUB;
use POSIX qw(log10);

class Genome::Model::Tools::Vcf::Convert::Snv::VarscanSomatic {
    is => 'Genome::Model::Tools::Vcf::Convert::Base',
    doc => 'Generate a VCF file from varscan somatic output'
};

sub help_synopsis {
    <<'HELP';
    Generate a VCF file from varscan somatic snv output
HELP
}

sub help_detail {
    <<'HELP';
    Parses the input file and creates a VCF containing all the snvs.
HELP
}

sub source {
    my $self = shift;
    return "VarscanSomatic";
}

sub parse_line {
    my $self = shift;
    my $line = shift;

    # TODO snv_qual is not used...
    my @col = split("\t",$line);
    my $chr = $col[0];
    my $pos = $col[1];
    my $ref = $col[2];
    my $variant = $col[3];
    my $normal_genotype = $col[7];
    my $tumor_genotype = $col[11];
    #total read depth
    my $normal_dp  = $col[4]+$col[5];
    my $tumor_dp = $col[8]+$col[9];
    #replace ambiguous/IUPAC bases with N in ref
    $ref =~ s/[^ACGTN\-]/N/g;

    # Get all of the possible ALT calls
    my @alt_alleles = Genome::Info::IUB->variant_alleles_for_iub($ref, $variant);
    push (@alt_alleles, Genome::Info::IUB->variant_alleles_for_iub($ref, $normal_genotype));
    push (@alt_alleles, Genome::Info::IUB->variant_alleles_for_iub($ref, $tumor_genotype));
    # Get unique alt alleles
    my %unique_alleles;
    @alt_alleles = grep(!$unique_alleles{$_}++, @alt_alleles);
    my $alt = join(",", @alt_alleles);

    #add the ref and alt alleles' positions in the allele array to the GT field
    my @tumor_alleles = Genome::Info::IUB->iub_to_alleles($tumor_genotype);
    my $tumor_gt = $self->generate_gt($ref, \@alt_alleles, \@tumor_alleles);

    #add the ref and alt alleles' positions in the allele array to the GT field
    my @normal_alleles = Genome::Info::IUB->iub_to_alleles($normal_genotype);
    my $normal_gt = $self->generate_gt($ref, \@alt_alleles, \@normal_alleles);

    # We do not have access to much of the normal information from somatic output
    #genotype quality (consensus quality)
    my $normal_gq = ".";
    my $tumor_gq = ".";
    # avg mapping quality ref/var
    my $normal_mq = ".";
    my $tumor_mq = ".";
    # avg mapping quality ref/var
    my $normal_bq = ".";
    my $tumor_bq = ".";
    # allele depth
    my $normal_ad =  ".";
    my $tumor_ad =  ".";
    # fraction of reads supporting alt
    $col[6] =~ s/\%// ;
    $col[10] =~ s/\%// ;
    my $normal_fa = $col[6]/100;
    my $tumor_fa = $col[10]/100;
    # vaq
    #edge case where a score of zero results in "inf"
    my $tumor_vaq;
    if($col[14] == 0){
        $tumor_vaq = 99;
    } else {
        $tumor_vaq = sprintf "%.0f", -10*log10($col[14]);
    }
    my $normal_vaq = ".";

    # Placeholder for later adjustment
    my $dbsnp_id = ".";
    my $qual = "."; # Can also be $tumor_vaq
    my $filter = "PASS";
    my $format = "GT:GQ:DP:BQ:MQ:AD:FA:VAQ";
    my $info = ".";
    my $tumor_sample_string = join (":", ($tumor_gt, $tumor_gq, $tumor_dp, $tumor_bq, $tumor_mq, $tumor_ad, $tumor_fa, $tumor_vaq));
    my $normal_sample_string = join (":", ($normal_gt, $normal_gq, $normal_dp, $normal_bq, $normal_mq, $normal_ad, $normal_fa, $normal_vaq));

    my $vcf_line = join("\t", $chr, $pos, $dbsnp_id, $ref, $alt, $qual, $filter, $info, $format, $normal_sample_string, $tumor_sample_string);

    return $vcf_line;
}

