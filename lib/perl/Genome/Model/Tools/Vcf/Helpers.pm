package Genome::Model::Tools::Vcf::Helpers;

use strict;
use warnings;

use above 'Genome';
use Genome::Info::IUB;

require Exporter;
our @ISA = qw/Exporter/;
our @EXPORT_OK = qw/
    convertIub
    genGT
/;


sub convertIub{
    my ($base) = @_;
    #deal with cases like "A/T" or "C/T"
    if ($base =~/\//){
        my @bases=split(/\//,$base);
        my %baseHash;
        foreach my $b (@bases){
            my $res = convertIub($b);
            my @bases2 = split(",",$res);
            foreach my $b2 (@bases2){
                $baseHash{$b2} = 0;
            }
        }
        return join(",",keys(%baseHash));
    }

    return join(',', Genome::Info::IUB->iub_to_bases($base)); #TODO: this is completely stupid.  make this not
}

# generate a GT line from a base and a list of all alleles at the position
sub genGT{
    my ($base, @alleles) = @_;
    my @bases = split(",",convertIub($base));
    if (@bases > 1){
        my @pos;
        push(@pos, (firstidx{ $_ eq $bases[0] } @alleles));
        push(@pos, (firstidx{ $_ eq $bases[1] } @alleles));
        return(join("/", sort(@pos)));
    } else { #only one base
        my @pos;
        push(@pos, (firstidx{ $_ eq $bases[0] } @alleles));
        push(@pos, (firstidx{ $_ eq $bases[0] } @alleles));
        return(join("/", sort(@pos)));
    }
}

1;
