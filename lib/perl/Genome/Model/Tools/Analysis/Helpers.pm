package Genome::Model::Tools::Analysis::Helpers;

use strict;
use warnings;

require Exporter;
our @ISA = qw(Exporter);
our @EXPORT_OK = qw(
    byBamOrder
    code_to_genotype
    code_to_genotype_returning_code
    commify
    flip_genotype
    is_heterozygous
    is_homozygous
    sort_genotype
);

sub byBamOrder ($$)
{
    my ($chrom_a, $pos_a) = split(/\t/, $_[0]);
    my ($chrom_b, $pos_b) = split(/\t/, $_[1]);

    $chrom_a =~ s/X/9\.1/;
    $chrom_a =~ s/Y/9\.2/;
    $chrom_a =~ s/MT/25/;
    $chrom_a =~ s/M/25/;
    $chrom_a =~ s/NT/99/;
    $chrom_a =~ s/[^0-9\.]//g;

    $chrom_b =~ s/X/9\.1/;
    $chrom_b =~ s/Y/9\.2/;
    $chrom_b =~ s/MT/25/;
    $chrom_b =~ s/M/25/;
    $chrom_b =~ s/NT/99/;
    $chrom_b =~ s/[^0-9\.]//g;

    $chrom_a <=> $chrom_b
    or
    $pos_a <=> $pos_b;
}

sub is_heterozygous
{
    my $gt = shift(@_);
    (my $a1, my $a2) = split(//, $gt);
    return(1) if($a1 ne $a2);
    return(0);
}

sub is_homozygous
{
    my $gt = shift(@_);
    (my $a1, my $a2) = split(//, $gt);
    return(1) if($a1 eq $a2);
    return(0);
}

sub flip_genotype
{
    my $gt = shift(@_);
    (my $a1, my $a2) = split(//, $gt);

    if($a1 eq "A")
    {
        $a1 = "T";
    }
    elsif($a1 eq "C")
    {
        $a1 = "G";
    }
    elsif($a1 eq "G")
    {
        $a1 = "C";
    }
    elsif($a1 eq "T")
    {
        $a1 = "A";
    }

    if($a2 eq "A")
    {
        $a2 = "T";
    }
    elsif($a2 eq "C")
    {
        $a2 = "G";
    }
    elsif($a2 eq "G")
    {
        $a2 = "C";
    }
    elsif($a2 eq "T")
    {
        $a2 = "A";
    }

    $gt = $a1 . $a2;
    $gt = sort_genotype($gt);
    return($gt);
}

sub sort_genotype
{
    my $gt = shift(@_);
    (my $a1, my $a2) = split(//, $gt);

    my @unsorted = ($a1, $a2);
    my @sorted = sort @unsorted;
    $a1 = $sorted[0];
    $a2 = $sorted[1];
    return($a1 . $a2);
}

sub code_to_genotype
{
    my $code = shift(@_);

    return("AA") if($code eq "A");
    return("CC") if($code eq "C");
    return("GG") if($code eq "G");
    return("TT") if($code eq "T");

    return("AC") if($code eq "M");
    return("AG") if($code eq "R");
    return("AT") if($code eq "W");
    return("CG") if($code eq "S");
    return("CT") if($code eq "Y");
    return("GT") if($code eq "K");

    return("NN");
}

sub code_to_genotype_returning_code {
    my $genotype = code_to_genotype(@_);
    return $genotype unless $genotype eq "NN";
    return shift;
}

sub commify
{
    local($_) = shift;
    1 while s/^(-?\d+)(\d{3})/$1,$2/;
    return $_;
}

1;
