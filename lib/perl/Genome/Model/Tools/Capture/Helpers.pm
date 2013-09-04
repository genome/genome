package Genome::Model::Tools::Capture::Helpers;

use strict;
use warnings;

use above 'Genome';

require Exporter;
our @ISA = qw(Exporter);
our @EXPORT_OK = qw(
    iupac_to_base
    fix_chrom
    sortByChrPos
    byChrPos
);

#############################################################
# IUPAC to base - convert IUPAC code to variant base
#
#############################################################

sub iupac_to_base
{
    (my $allele1, my $allele2) = @_;

    return($allele2) if($allele2 eq "A" || $allele2 eq "C" || $allele2 eq "G" || $allele2 eq "T");

    if($allele2 eq "M")
    {
        return("C") if($allele1 eq "A");
        return("A") if($allele1 eq "C");
    }
    elsif($allele2 eq "R")
    {
        return("G") if($allele1 eq "A");
        return("A") if($allele1 eq "G");
    }
    elsif($allele2 eq "W")
    {
        return("T") if($allele1 eq "A");
        return("A") if($allele1 eq "T");
    }
    elsif($allele2 eq "S")
    {
        return("C") if($allele1 eq "G");
        return("G") if($allele1 eq "C");
    }
    elsif($allele2 eq "Y")
    {
        return("C") if($allele1 eq "T");
        return("T") if($allele1 eq "C");
    }
    elsif($allele2 eq "K")
    {
        return("G") if($allele1 eq "T");
        return("T") if($allele1 eq "G");
    }

    return($allele2);
}

sub fix_chrom
{
    my $chrom = shift;
    $chrom =~ s/chr// if(substr($chrom, 0, 3) eq "chr");
    $chrom =~ s/[^0-9XYMNT\_random]//g;

    return($chrom);
}


sub byChrPos
{
    (my $chrom_a, my $pos_a) = split(/\t/, $a);
    (my $chrom_b, my $pos_b) = split(/\t/, $b);

    $chrom_a =~ s/X/23/;
    $chrom_a =~ s/Y/24/;
    $chrom_a =~ s/MT/25/;
    $chrom_a =~ s/M/25/;
    $chrom_a =~ s/[^0-9]//g;

    $chrom_b =~ s/X/23/;
    $chrom_b =~ s/Y/24/;
    $chrom_b =~ s/MT/25/;
    $chrom_b =~ s/M/25/;
    $chrom_b =~ s/[^0-9]//g;

    $chrom_a <=> $chrom_b
    or
    $pos_a <=> $pos_b;
}

sub sortByChrPos{
    return sort byChrPos @_;
}

1;
