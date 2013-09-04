#!/usr/bin/env genome-perl

use warnings;
use strict;
use IO::File;
use above 'Genome';
use Genome::Info::IUB;

# arg 0 = snps (5 col (chr, pos, pos, ref, var)
# arg 1 = readcounts


my %refHash;
my %varHash;


#read in all the snvs and hash both the ref and var allele by position
my $inFh = IO::File->new( $ARGV[0] ) || die "can't open file\n";
while( my $line = $inFh->getline )
{
    chomp($line);
    my @fields = split("\t",$line);
    $refHash{$fields[0] . "|" . $fields[1]} = $fields[3];
    $varHash{$fields[0] . "|" . $fields[1]} = $fields[4]
}


#convert weird bases to lists
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
    };

#
sub matchIub{
    my ($allele,$ref,$var) = @_;
    my @variubs = split(",",convertIub($var));
    my @refiubs = split(",",convertIub($ref));
    foreach my $i (@variubs){
	unless (grep {$_ eq $i} @refiubs) {
	    if ($allele eq $i){
		return 1;
	    }
	}
    }
    return 0;
}


#read in the bam-readcount file
my $inFh2 = IO::File->new( $ARGV[1] ) || die "can't open file\n";
while( my $line = $inFh2->getline )
{
    chomp($line);
    my ($chr, $pos, $ref, $depth, @counts) = split("\t",$line);

    my $ref_count = 0;
    my $var_count = 0;
    my $knownRef;
    my $knownVar;
    
    #for each base at that pos
    foreach my $count_stats (@counts) {
	my ($allele, $count, $mq, $bq) = split /:/, $count_stats;

	# skip if it's not in our list of snvs
	next unless (exists($refHash{$chr . "|" . $pos}) && exists($varHash{$chr . "|" . $pos}));
	
	#look up the snv calls at this position
	$knownRef = $refHash{$chr . "|" . $pos};
	$knownVar = $varHash{$chr . "|" . $pos};

	# assume that the ref call is ACTG, not iub 
	# (assumption looks valid in my files)
	if ($allele eq $knownRef){
	    $ref_count += $count;
	}
	
	# if this base is included in the IUB code for
	# for the variant, (but doesn't match the ref)
	if (matchIub($allele,$knownRef,$knownVar)){
	    $var_count += $count;
	}
    }

    my $var_freq = 0;
    if ($depth ne '0') {
        $var_freq = $var_count/$depth * 100;
    }

    #output
    print "$chr\t$pos\t$knownRef\t$knownVar\t$ref_count\t$var_count\t";
    printf("%.2f",$var_freq);
    print "\n";
}
