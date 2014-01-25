#!/usr/bin/env genome-perl
# This script will read Cufflink's GTF and get the junctions used in the transcripts
# Rodrigo Goya
use strict;
use warnings;

my $GTF_INFILE = $ARGV[0] || "-";
my $line_number = 0;
my %trs;
open(GTF,"<$GTF_INFILE") || die("ERROR: could not read GTF fileL '$GTF_INFILE'\n");
while(<GTF>) {
  chomp; s///;
  if(/^#/) { next; }
  #chr1	Cufflinks	transcript	898801	899131	1000	+	.	gene_id "CUFF.95"; transcript_id "CUFF.95.1"; FPKM "10.1796906583"; frac "1.000000"; conf_lo "5.934466"; conf_hi "14.424915"; cov "8.935252";
  #chr1	Cufflinks	exon	898801	898883	1000	+	.	gene_id "CUFF.95"; transcript_id "CUFF.95.1"; exon_number "1"; FPKM "10.1796906583"; frac "1.000000"; conf_lo "5.934466"; conf_hi "14.424915"; cov "8.935252";
  #chr1	Cufflinks	exon	899076	899131	1000	+	.	gene_id "CUFF.95"; transcript_id "CUFF.95.1"; exon_number "2"; FPKM "10.1796906583"; frac "1.000000"; conf_lo "5.934466"; conf_hi "14.424915"; cov "8.935252";
  my @f = split(/\t/, $_);
  my %flags;
  foreach my $flag (split(/;\s*/, $f[$#f])) {
    $flag =~ m/^\s*(\S+)\s*"(\S+)"/;
    $flags{$1}= $2;
  }
  $line_number++;
  if($f[2] eq "transcript"){
    if(exists($trs{$flags{transcript_id}})) {
      die("ERROR: re-definition of cufflinks transcript_id = $flags{transcript_id},  skipping\n");
    }
    $trs{$flags{transcript_id}}{coords}{chr} = $f[0];
    $trs{$flags{transcript_id}}{coords}{start} = $f[3];
    $trs{$flags{transcript_id}}{coords}{end} = $f[4];
    $trs{$flags{transcript_id}}{coords}{strand} = $f[6];
    $trs{$flags{transcript_id}}{exon_num} = 0;
    $trs{$flags{transcript_id}}{fpkm} = $flags{FPKM};
    #print STDERR $line_number." transcript $flags{transcript_id}\n";
  }elsif($f[2] eq "exon"){
    if(!exists($trs{$flags{transcript_id}})) {
      die("ERROR: transcript $flags{transcript_id} has not yet been defined\n");
  }
  $trs{$flags{transcript_id}}{exons}[$flags{exon_number}-1]{chr} = $f[0];
  $trs{$flags{transcript_id}}{exons}[$flags{exon_number}-1]{start} = $f[3];
  $trs{$flags{transcript_id}}{exons}[$flags{exon_number}-1]{end} = $f[4];
  $trs{$flags{transcript_id}}{exon_num}++;
  #print STDERR $line_number." exon $flags{exon_number}\n";
  }else{
    warn("WARNING: unknown feature in line $line_number\n");
    next;
  }
}
close(GTF);

my %stats;
my %juncs;
foreach my $tr (keys %trs){
  $stats{exon_num}[$trs{$tr}{exon_num}]++;
  if($trs{$tr}{exon_num} > 1) {
    my $left = -1;
    for(my $i = 0; $i < $trs{$tr}{exon_num}; $i++) {
      if($left != -1) {
        print "$trs{$tr}{coords}{chr}:$left-$trs{$tr}{exons}[$i]{start}\t$trs{$tr}{fpkm}\n";
      } 
      $left = $trs{$tr}{exons}[$i]{end};
    }
  }
}

for(my $i = 0; $i < scalar(@{$stats{exon_num}}); $i++) {
  print STDERR "$i exon transcripts: $stats{exon_num}[$i]\n";
}

exit();


