#!/usr/bin/env perl
#Written by Malachi Griffith
#The purpose of this script is to get data for a chromosome segment specified by coordinates

#Make sure you select a version of the ensembl API that matches the database and genome build you wish to use!!
#For example to use homo_sapiens_core_35_35h  (Ensembl version 35, genome version 35h which equals hg17)

use strict;
use warnings;
use Data::Dumper;
use Getopt::Long;
use Term::ANSIColor qw(:constants);

#Initialize command line options
my $ensembl_api_version = ''; #Version of EnsEMBL to use
my $ensembl_host = '';
my $ensembl_user = '';
my $species = '';
my $name1 = '';
my $chr1 = '';
my $start1 = '';
my $end1 = '';
my $name2 = '';
my $chr2 = '';
my $start2 = '';
my $end2 = '';
my $outdir = '';

GetOptions ('ensembl_api_version=s'=>\$ensembl_api_version,
	          'ensembl_host=s'=>\$ensembl_host, 'ensembl_user=s'=>\$ensembl_user,
	          'species=s'=>\$species,
            'chr1=s'=>\$chr1, 'start1=i'=>\$start1, 'end1=i'=>\$end1, 'name1=s'=>\$name1, 
            'chr2=s'=>\$chr2, 'start2=i'=>\$start2, 'end2=i'=>\$end2, 'name2=s'=>\$name2,
            'outdir=s'=>\$outdir);

#Provide instruction to the user
unless (($ensembl_api_version =~ /^\d+/) && $ensembl_host && $ensembl_user && $species && $chr1 && $start1 && $end1 && $name1 && $chr2 && $start2 && $end2 && $name2 && $outdir){
  print GREEN, "\n\nExample: getSvBreakpointDnaSequence.pl  --ensembl_api_version=70  --ensembl_host='ensembldb.ensembl.org'  --ensembl_user='anonymous'  --species='Human'  --chr1=22  --start1=41536000  --end1=41536050  --name1='EP300_5prime'  --chr2=12  --start2=6785285  --end2=6785235  --name2='ZNF384_3prime' --outdir=/gscmnt/sata132/techd/mgriffit/all1/sv_fusion_primer_design/\n", RESET;
  print GREEN, "\n\nTo get the reverse complement of one or both sequences, make --start larger than --end\n", RESET;
  print RED, "\nBasic option(s) missing or incorrect format\n\n", RESET;
  exit();
}

if ($ensembl_api_version =~ /^\d+/){
  if($ensembl_api_version eq "70"){
    unshift(@INC, "/gscmnt/sata206/techd/ensembl_api/v70/ensembl/modules");
  }else{
    print RED, "\nEnsEMBL API version: $ensembl_api_version is not defined, modify script before proceeding\n\n", RESET;
    exit();
  }
}else{
  print RED, "\nEnsEMBL API version format: $ensembl_api_version not understood!\n\n", RESET;
  exit();
}
use lib "/gscmnt/sata206/techd/bioperl-live";    #Bioperl
require Bio::EnsEMBL::Registry;  #Use for remote connections over the web


#Check the outdir
$outdir .= "/" unless ($outdir =~ /\/$/);
die "\n\nNot a valid outdir: $outdir\n\n" unless (-e $outdir && -d $outdir);

# First, we need to import all Perl modules that we will be using. 
#Since we need a connection to an Ensembl database through the Registry we first have to import the Registry module which we use to establish this connection.
my $registry = 'Bio::EnsEMBL::Registry';
$registry->load_registry_from_db(
    -host => $ensembl_host, # alternatively 'useastdb.ensembl.org'
    -user => $ensembl_user
);

my $slice_adaptor = $registry->get_adaptor($species, 'Core', 'Slice');

#Get the first sequence
my $slice1;
if ($start1 < $end1){
  $slice1 = $slice_adaptor->fetch_by_region('chromosome', $chr1, $start1, $end1);
}else{
  $slice1 = $slice_adaptor->fetch_by_region('chromosome', $chr1, $end1, $start1);
}
my $coord_sys1 = $slice1->coord_system()->name();
my $seq_region1 = $slice1->seq_region_name();
my $seq_start1 = $slice1->start();
my $seq_end1 = $slice1->end();
my $strand1 = $slice1->strand();
my $sequence1 = $slice1->seq();
my $sequence1_rc = reverse($sequence1);
$sequence1_rc =~ tr/ACGTacgt/TGCAtgca/;
print BLUE, "\n\nSlice: coord_sys1=$coord_sys1\tseq_region1=$seq_region1\tseq_start1=$seq_start1\tseq_end1=$seq_end1\tstrand1=$strand1\n$sequence1\n\n", RESET;


#Get the second sequence
my $slice2;
if ($start2 < $end2){
  $slice2 = $slice_adaptor->fetch_by_region('chromosome', $chr2, $start2, $end2);
}else{
  $slice2 = $slice_adaptor->fetch_by_region('chromosome', $chr2, $end2, $start2);
}
my $coord_sys2 = $slice2->coord_system()->name();
my $seq_region2 = $slice2->seq_region_name();
my $seq_start2 = $slice2->start();
my $seq_end2 = $slice2->end();
my $strand2 = $slice2->strand();
my $sequence2 = $slice2->seq();
my $sequence2_rc = reverse($sequence2);
$sequence2_rc =~ tr/ACGTacgt/TGCAtgca/;
print BLUE, "\n\nSlice: coord_sys2=$coord_sys2\tseq_region2=$seq_region2\tseq_start2=$seq_start2\tseq_end2=$seq_end2\tstrand2=$strand2\n$sequence2", RESET;

#Piece the two sequences together and separate by [N] to ensure that primer design by PRIMER3 will flank the break-point
my $fusion_seq = '';
my $fusion_name = '';
if ($start1 < $end1){
  $fusion_seq .= $sequence1;
  $fusion_name .= $name1;
}else{
  $fusion_seq .= $sequence1_rc;
  $fusion_name .= $name1 . "(rc)";
}
$fusion_seq .= "[N]";
$fusion_name .= "-";
if ($start2 < $end2){
  $fusion_seq .= $sequence2;
  $fusion_name .= $name2;
}else{
  $fusion_seq .= $sequence2_rc;
  $fusion_name .= $name2 . "(rc)";
}

my $fusion_file_name1 = $name1 . "-" . $name2 . ".fa";
my $fusion_outfile1 = $outdir . $fusion_file_name1;

#Print out the complete sequence
open (OUT1, ">$fusion_outfile1") || die "\n\nCould not open output file for writing: $fusion_outfile1\n\n";
my $outstring = ">$fusion_name\n$fusion_seq", RESET;
#print BLUE, "\n\n$outstring", RESET;
print OUT1 "$outstring";
close (OUT1);

#Print out a truncated version that has less flank on either side of the fusion
my $flank = 1000;
my $fusion_file_name2 = $name1 . "-" . $name2 . "." . $flank ."bases.fa";
my $fusion_outfile2 = $outdir . $fusion_file_name2;
if ($fusion_seq =~ /(\w{$flank}\[N\]\w{$flank})/){
  my $fusion_seq_flank = $1;
  open (OUT2, ">$fusion_outfile2") || die "\n\nCould not open output file for writing: $fusion_outfile2\n\n";
  my $outstring = ">$fusion_name\n$fusion_seq_flank", RESET;
  #print BLUE, "\n\n$outstring", RESET;
  print OUT2 "$outstring";
  close (OUT2);
}else{
  print "\n\nWARNING: Not enough sequence to produce the shorter version with only $flank bases of flank on each side";
}
print BLUE, "\n\nPrinted SV primer design files: \n$fusion_outfile1\n$fusion_outfile2\n\n", RESET;

exit();

