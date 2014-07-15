#!/usr/bin/env genome-perl

# AUTHOR: Joseph Fass
# LAST REVISED: January 2010
# 
# The Bioinformatics Core at UC Davis Genome Center
# http://bioinformatics.ucdavis.edu
# Copyright (c) 2010 The Regents of University of California, Davis Campus.
# All rights reserved.

use strict;
use Getopt::Std;

my $usage = "\nusage: $0 -f <unmapped bam file> [-q # -o # -l #]\n\n".
    "Trims sequences from BAM file a la Heng Li's bwa '-q' option\n\n".
    "\t-f    \tUnmapped bam file to trim\n".
    "\t-q int\tQuality threshold (default:  20).  To trim q2, then set -q to 3\n".
    "\t-o int\tFastq offset value (default:  33).  Sanger fastq is default, change to 65 for Illumina.\n".
    "\t-l int\tLength cutoff of trimmed reads (default:  60).\n".
    "\t      \tIf both pairs are > length cutoff, then they will be put into read1/read2 fastq.\n".
    "\t      \tIf only 1 read is > length cutoff, then it will be put into a singleton file.\n\n".
    "For specifics about quality trimming bwa style (see bwa's bwa_trim_read() function in bwaseqio.c)\n\n";

our($opt_q,$opt_o,$opt_f,$opt_l,$opt_h);
getopts('hq:o:f:l:');
if (!defined($opt_q) or !($opt_q =~ m/^[0-9]+$/)) {$opt_q = 20;}
if (!defined($opt_o) or !($opt_o =~ m/^[0-9]+$/)) {$opt_o = 33;}
if (!defined($opt_l) or !($opt_l =~ m/^[0-9]+$/)) {$opt_l = 60;}

if (!$opt_f || $opt_h) {die $usage;}

chomp(my $samtools = `which samtools`);
if (!$samtools) {
    die "Samtools not found.  Please check path.\n";
}

my $pos;  my $maxPos;  my $area;  my $maxArea;

my $read1_file = $opt_f;
my $read2_file = $opt_f;
my $singletonFile = $opt_f;
$read1_file =~ s/bam/trimmed.1.fastq/;
$read2_file =~ s/bam/trimmed.2.fastq/;
$singletonFile =~ s/bam/trimmed.singleton.fastq/;
my $quality_reads = 0;
my $quality_bases = 0;

open (FILE, "$samtools view $opt_f |") or die "Can't open bam for reading, $opt_f:  $!\n";
open (READ1, ">$read1_file") or die "Can't open read1 file, $read1_file\n";
open (READ2, ">$read2_file") or die "Can't open read2 file, $read2_file\n";
open (SINGLETON, ">$singletonFile") or die "Can't open singleton file, $singletonFile:  $!\n";

while (<FILE>) {
    chomp;
    my $read1 = $_;
    my @array = split/\t/, $read1;
    my $read1_name = $array[0]."/1";
    my $read1_s = $array[9];
    my $read1_q = $array[10];
    my $read1_length = &checkPos($read1_q);
    my $read1_seq = substr($read1_s,0,$read1_length);  
    my $read1_qual = substr($read1_q,0,$read1_length);
    my $read2 = <FILE>;
    undef(@array);
    @array = split/\t/, $read2;
    my $read2_name = $array[0]."/2";
    my $read2_s = $array[9];
    my $read2_q = $array[10];
    my $read2_length = &checkPos($read2_q);
    my $read2_seq = substr($read2_s,0,$read2_length);  
    my $read2_qual = substr($read2_q,0,$read2_length);
    if ($read1_length >= $opt_l && $read2_length >= $opt_l) {
	print READ1 "\@$read1_name\n$read1_seq\n\+$read1_name\n$read1_qual\n";
	print READ2 "\@$read2_name\n$read2_seq\n\+$read2_name\n$read2_qual\n";
	$quality_reads += 2;
	$quality_bases += $read1_length;
	$quality_bases += $read2_length;
    } elsif ($read1_length < $opt_l && $read2_length >= $opt_l) {
	print SINGLETON "\@$read2_name\n$read2_seq\n\+$read2_name\n$read2_qual\n";
	$quality_reads++;
	$quality_bases += $read2_length;
    } elsif ($read2_length < $opt_l && $read1_length >= $opt_l) {
	print SINGLETON "\@$read1_name\n$read1_seq\n\+$read1_name\n$read1_qual\n";
	$quality_reads++;
	$quality_bases += $read1_length;
    }

}

print "Quality reads: $quality_reads\n";
print "Quality bases: $quality_bases\n";

close FILE;
close READ1;
close READ2;
close SINGLETON;

sub checkPos {
    my $q = shift;
    $pos = length($q);
    $maxPos = $pos;
    $area = 0;
    $maxArea = 0;
    while ($pos>0 && $area>=0) {
	$area += $opt_q - (ord(substr($q,$pos-1,1))-$opt_o);
	if ($area > $maxArea) {
	    $maxArea = $area;
	    $maxPos = $pos;
	}
	$pos--;
	
    }  
    $maxPos = ($pos == 0 ? 1 : $maxPos);
    return $maxPos;
}

