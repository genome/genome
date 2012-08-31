#!/usr/bin/env genome-perl -w
# -*- Perl -*-




#
#
#  processSvReadRemapOutFiles.pl
#
#
# Parse output from 'svCaptureValidation.pl'
#
# Look at the total number of reads from a given region
#             total number of SV-supporting reads
#
# There must be exactly one normal sample and can be any number of other (tumor) samples
#
# Print out all lines with p-value of tumor to normal SV read comparison
# Consider SV to be somatic if:
#        any tumor sample has >= $minTumorSvReads and >= $minTotalReads and
#        normal sample has <= $maxNormalSvReads and
#        a chi-square test says number of that tumor SV reads is significantly different than 
#        number of normal SV reads given the coverage (i.e. total number of reads from region)
#
#  


use strict;
use Carp;

my $minTumorSvReads = 5;
my $maxNormalSvReads = 10;
my $minTotalReads = 10;
my $maxPvalue = 0.05;

my ( $svFile, $line, @entireFile, $normalTotal, $normalSv, $tumorTotal, $tumorSv, $pvalue, 
    %totalReads, %svReads, $hasMinTumorSvReads, $somatic, $sampleName, @cols, $entry, $filtered, 
    );


$svFile = $ARGV[0];
(-e $svFile) || die "Annotated SV file not found: '$svFile'.  Provide full path to file";


open(IN, "< $svFile") || die "Could not open '$svFile': $!";
@entireFile = <IN>;
close IN;
foreach $line (@entireFile) {
    chomp $line;
    print "$line";
    if ( $line =~ /\#/ ) { print "\n"; next; }

    %totalReads = %svReads = ();
    $normalTotal = undef;
    $hasMinTumorSvReads = 0;
    @cols = split /\s+/, $line;
    foreach $entry ( @cols ) {
	if ( $entry =~ /normal\.totalReads\:(\d+)/i ) {
	    (!defined $normalTotal) || die "More than one normal sample in '$line'";
	    $normalTotal = $1;
	} elsif ( $entry =~ /normal\.\w*svReadCount\:(\d+)/i ) {
	    $normalSv = $1;
	} elsif ( $entry =~ /(\S+)\.totalReads\:(\d+)/ ) {
	    $totalReads{$1} = $2;
	} elsif ( $entry =~ /(\S+)\.\w*svReadCount\:(\d+)/i ) {
	    $svReads{$1} = $2;
	    if ( $2 >= $minTumorSvReads ) { $hasMinTumorSvReads = 1; }
	}
    }

    # Calculate p-values for all samples vs. normal
    $somatic = $filtered = "-";
    foreach $sampleName (keys %totalReads) {
	$tumorTotal = $totalReads{$sampleName};
	$tumorSv = $svReads{$sampleName};
	(defined $tumorSv) || die "Did not get number of SV reads for $sampleName in '$line'";
	if ( $normalSv == 0 && $tumorSv == 0 ) {
	    # pvalue is undefined
	    print "\t$sampleName"."_normal:undef";
	    next;
	}	
	$pvalue = pValue($normalTotal, $normalSv, $tumorTotal, $tumorSv);
	print "\t$sampleName"."_normal:$pvalue";
	(defined $pvalue) || die "p-value is undefined for ($normalTotal, $normalSv, $tumorTotal, $tumorSv) in '$line'";
	if ( $pvalue <= $maxPvalue && $tumorSv > $normalSv ) { $somatic = "somatic"; }
    }

    # Calculate p-values comparing all tumor samples to all other tumor samples
    if ( scalar(keys %totalReads) >= 2 ) {
	my @tumorNames = keys %totalReads;
	for ( my $i = 0; $i < $#tumorNames; $i++ ) {
	    my $j = $i + 1;
	    my $firstSample = $tumorNames[$i];
	    my $secondSample = $tumorNames[$j];
	    if ( $svReads{$firstSample} == 0 && $svReads{$secondSample} == 0 ) {
		print "\t$firstSample"."_"."$secondSample:undef";
		next;
	    }
	    $pvalue = pValue($totalReads{$firstSample}, $svReads{$firstSample}, $totalReads{$secondSample}, $svReads{$secondSample});
	    print "\t$firstSample"."_"."$secondSample:$pvalue";
	}
    }
		
	    

    # Filter out if there are not enough total normal reads or too many normal SV reads
    # or none of the 'tumor' samples have the minimum number of SV-supporting reads
    if ( $somatic eq "somatic" && 
	($normalTotal < $minTotalReads || $normalSv > $maxNormalSvReads || !$hasMinTumorSvReads)
	) { $filtered = "filtered"; }

    # Print out if any of the samples are somatic when compared to normal
    print "\t$somatic\t$filtered\n";

}

sub pValue {
    # The world's most inefficient way to get a pvalue
    # Input: total normal reads, normal SV reads, total tumor reads, tumor SV reads
    # Return: pvalue from Fisher's Exact test to see if the SV reads from
    #         tumor and normal are significantly different given the total number of
    #         reads from each

    my ($normalTotal, $normalSv, $tumorTotal, $tumorSv) = @_;
    
    my $pvalue;
    my $rFile = "/tmp/rFile".rand();
    open(OUT, "> $rFile") || die "Could not open '$rFile': $!";
    print OUT "x=matrix(c($tumorTotal,$tumorSv,$normalTotal,$normalSv),ncol=2)\nfisher.test(x)\n";
    close OUT;
    open(ROUT, "R --silent --slave <  $rFile |");
    my @rOutput = <ROUT>;
    close ROUT;
    foreach (@rOutput) { 
	if ( $_ =~ /p-value\s+[=<]\s+(\S+)/ ) { $pvalue = $1; }
    }
    unlink $rFile;

    return $pvalue;
}


sub BAK_pValue {
    # The world's most inefficient way to get a pvalue
    # Input: total normal reads, normal SV reads, total tumor reads, tumor SV reads
    # Return: pvalue from chi-square test to see if the SV reads from
    #         tumor and normal are significantly different given the total number of
    #         reads from each

    my ($normalTotal, $normalSv, $tumorTotal, $tumorSv) = @_;
    
    my $pvalue;
    my $rFile = "/tmp/rFile".rand();
    open(OUT, "> $rFile") || die "Could not open '$rFile': $!";
    print OUT "x=matrix(c($tumorTotal,$tumorSv,$normalTotal,$normalSv),ncol=2)\nchisq.test(x)\n";
    close OUT;
    open(ROUT, "R --silent --slave <  $rFile |");
    my @rOutput = <ROUT>;
    close ROUT;
    foreach (@rOutput) { 
	if ( $_ =~ /p-value\s+[=<]\s+(\S+)/ ) { $pvalue = $1; }
    }
    unlink $rFile;

    return $pvalue;
}
