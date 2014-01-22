#!/usr/bin/env genome-perl
#Written by Nate Dees and modified by Malachi Griffith

use warnings;
use strict;
use IO::File;
use above 'Genome';
use Getopt::Long;
use Term::ANSIColor qw(:constants);
use Data::Dumper;

my $script_dir;
use Cwd 'abs_path';
if (abs_path($0) =~ /(.*\/).*\.pl/){
  $script_dir = $1;
}

my $sites_file = '';
my $bam_list = '';
my $reference_fasta = '';
my $output_file = '';

GetOptions ('sites_file=s'=>\$sites_file, 'bam_list=s'=>\$bam_list, 'reference_fasta=s'=>\$reference_fasta, 'output_file=s'=>\$output_file);

my $usage=<<INFO;

  Example usage: 
  
  give_me_readcounts.pl  --sites_file=allsnvs.hq.novel.tier123.v2.bed.adapted  --bam_list="Tumor:/gscmnt/ams1183/info/build_merged_alignments/merged-alignment-blade14-2-15.gsc.wustl.edu-mcordes-4636-111286100/111286100.bam,Normal:/gscmnt/ams1181/info/build_merged_alignments/merged-alignment-blade13-4-13.gsc.wustl.edu-mcordes-31373-111136633/111136633.bam"  --reference_fasta=/gscmnt/sata420/info/model_data/2857786885/build102671028/all_sequences.fa  --output_file=allsnvs.hq.novel.tier123.v2.bed.adapted.readcounts
  
  Intro:
  This script gets read counts from a series of BAMs for the positions specified in the sites file

  Details:
  --sites_file                    PATH.   List of chromosome sites to get BAM read counts for. Format: (1	16949959	16949959	G	R)
  --bam_list                      STRING. List of BAM files to interrogate at these positions
  --reference_fasta               PATH.   The genome reference fasta file that was used to create the BAM 
  --output_file                   PATH    The file to write results to

INFO

unless ($sites_file && $bam_list && $reference_fasta && $output_file){
  print GREEN, "$usage", RESET;
  exit();
}

#if(@ARGV<1){die "$0  adapted_sites_file  bam_file(s)(type:file,type:file...)  build(36or37)  output_file\n";}
#my ($sites,$bams,$build,$output) = @ARGV;
my $sites = $sites_file;
my $bams = $bam_list;
my $output = $output_file;
my @bams = split /,/,$bams;
my @types;
my $sitesfh = new IO::File $sites,"r";
my $outfh = new IO::File $output,"w";

#headers
#print $outfh "#",$bams,"\n";
print $outfh "#Chr\tStart\tStop";
for my $bam (@bams) {
    my ($type,$bam) = split /:/,$bam;
    push @types,$type;
    print $outfh "\t$type" . "_Ref\t$type" . "_Var\t$type" . "_Ref_Count\t$type" . "_Var_Count\t$type" . "_Var_Freq";
}
print $outfh "\n";

my $site_list_for_readcount = Genome::Sys->create_temp_file_path();
print `cut -f 1-3 $sites > $site_list_for_readcount`;

#begin to load output hash
my %output;
while (my $line = $sitesfh->getline) {
    chomp $line;
    my ($chr,$start,$stop) = split /\t/,$line;
    $output{$chr}{$start} = {};
}
$sitesfh->close;

#cycle through bams to get readcounts
for my $bam (@bams) {
  my $temp_readcount_file = Genome::Sys->create_temp_file_path();
  my ($type,$bam) = split /:/,$bam;
  
  my $rv = Genome::Model::Tools::Sam::Readcount->execute(
    bam_file => $bam,
    minimum_mapping_quality => 1,
    output_file => $temp_readcount_file,
    reference_fasta => $reference_fasta,
    region_list => $site_list_for_readcount,
  );
  #my $cmd = "bam-readcount -q 1 -l $site_list_for_readcount -f $reference_fasta $bam > $temp_readcount_file";
  #Genome::Sys->shellcmd(cmd => $cmd);

  my $statscmd = "$script_dir/miller-bam-readcount-to-stats.noheader.pl $sites $temp_readcount_file |";
  open(STATS,$statscmd) or die "Couldn't open stats command: $!";
  while (my $line = <STATS>) {
    chomp $line;
    my ($chr,$pos,$stats) = split /\t/,$line,3;
    $output{$chr}{$pos}{$type} = $stats;
  }
}

#print output
for my $chr (sort keys %output) {
  for my $pos (sort keys %{$output{$chr}}) {
     print $outfh "$chr\t$pos\t$pos";
     for my $type (@types) {
       if (exists $output{$chr}{$pos}{$type}) {
         print $outfh "\t".$output{$chr}{$pos}{$type};
       }else{
         print $outfh "\tNA\tNA\tNA\tNA\tNA";     
       }
     }
     print $outfh "\n";  
  }
}

exit();


