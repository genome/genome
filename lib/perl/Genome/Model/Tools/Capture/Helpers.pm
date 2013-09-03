package Genome::Model::Tools::Capture::Helpers;

use strict;
use warnings;

use above 'Genome';

require Exporter;
our @ISA = qw(Exporter);
our @EXPORT_OK = qw(iupac_to_base);

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

1;
