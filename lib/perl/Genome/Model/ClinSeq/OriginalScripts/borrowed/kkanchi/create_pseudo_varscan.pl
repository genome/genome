#!/usr/bin/env genome-perl

use warnings;
use strict;
use above 'Genome';
use IO::File;

#arg 0 = snp file
#arg 1 = both tumor and normal readcounts (formatted)

if(@ARGV<1){die "$0\tsnps_file\ttumor_normal_readcounts\n";}

my ($snps,$readcounts) = @ARGV;

my %tumHash;
#my %normHash;
#read tumor
my $rcFh = IO::File->new( $readcounts ) || die "can't open file\n";
while( my $line = $rcFh->getline )
{
    chomp($line);
    my @fields = split("\t",$line);
    $tumHash{$fields[0] . "|" . $fields[1]} = $line;
}
$rcFh->close;

my $snpsFh = IO::File->new( $snps ) || die "can't open file\n";
while( my $line = $snpsFh->getline )
{
    chomp($line);
    my @fields = split("\t",$line);
    if(exists($tumHash{$fields[0] . "|" . $fields[1]})){
	#if(exists($normHash{$fields[0] . "|" . $fields[1]})){
	    my @tum = split("\t",$tumHash{$fields[0] . "|" . $fields[1]});
	    #my @norm = split("\t",$normHash{$fields[0] . "|" . $fields[1]});

	    print join("\t",(@fields[0..1],@fields[3..4])) . "\t";
	    #print join("\t",@norm[5..7]) . "\t";
	    #print "NULL\t";
	    print join("\t",@tum[10..12]) . "\t";
            print "NULL\t";
            print join("\t",@tum[5..7]) . "\t";
	    print "NULL\tSomatic\tNULL\tNULL\tNULL\tNULL\tNULL\tNULL\n";	   
	
    }
}
