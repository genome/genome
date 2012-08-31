#!/usr/bin/env genome-perl
# Written by Malachi Griffith

# Determines the following for all observed junctionsn a file 
#1.) Junctions that can be anchored to a known junction (matching donor OR acceptor) - consider strand here 
#    - If the Donor matches to a gene on the positive strand, the non matching Acceptor should have a higher chromosome coordinate
#    - If the Donor matches to a gene on the negative strand, the non matching Acceptor should have a lower chromosome coordinate
#    - If the Acceptor matches to a gene on the positive strand, the non-matching Donor should have a lower chromosome coordinate
#    - If the Acceptor matches to a gene on the negative strand, the non-matching Donor should have a higher chromosome coordinate
#2.) Junctions that could not be anchored at all.  
#3.) Exon skipping 
#    - Number of known splice sites skipped
#    - Number of known exon clusters skipped
#4.) Transcript IDs of matching transcripts (where both Acceptor and Donor match, or if no match, then the corresponding anchored transcripts)...
#    - CCDS, Ensembl, MGC, Refseq, UCSC, Vega
#5.) Gene IDs of matching transcripts (non-redundant list)
#6.) Determine the number of exons skipped and exon donors/acceptors skipped by each junction

use strict;
use warnings;
use Data::Dumper;
use Getopt::Long;
use Term::ANSIColor qw(:constants);
use File::Basename;
use IO::File;
use Benchmark;
use above 'Genome'; #remove the 'above' when we turn this into a module
use Genome::Model::ClinSeq::Util qw(:all);

my $gene_annotation_dir = '';
my $gene_annotation_set = '';
my $gene_name_type = '';
my $obs_junction_file = '';
my $bedtools_bin_dir = '';
my $working_dir = '';
my $verbose = 0;

GetOptions ('obs_junction_file=s'=>\$obs_junction_file, 'verbose=i'=>\$verbose, 
            'bedtools_bin_dir=s'=>\$bedtools_bin_dir, 'working_dir=s'=>\$working_dir,
            'gene_annotation_dir=s'=>\$gene_annotation_dir, 'gene_annotation_set=s'=>\$gene_annotation_set, 'gene_name_type=s'=>\$gene_name_type);

my $usage=<<INFO;

  Example usage: 
  
  annotateObservedJunctions.pl  --obs_junction_file='/gscmnt/sata132/techd/mgriffit/hgs/hg1/qc/tophat/junctions.strand.junc'  --bedtools_bin_dir='/gsc/bin/'  --working_dir='/gscmnt/sata132/techd/mgriffit/hgs/hg1/qc/tophat/'  --gene_annotation_dir='/gscmnt/sata132/techd/mgriffit/reference_annotations/hg19/'  --gene_annotation_set='ALL'  --gene_name_type='symbol'

  --gene_annotation_dir   Specify the path to gene annotation files
  --gene_annotation_set   Specify name of the gene annotation set within this dir
  --gene_name_type        Specify the type of gene names to be used (determines with reference junction file will be used) using: --gene_name_type=symbol OR --gene_name_type=id
  --bedtools_bin_dir      Specify the path to the BEDTools binary dir
  --working_dir           Specify the base analysis dir.
  --verbose               To display more output, set to 1

INFO


unless ($obs_junction_file && $bedtools_bin_dir && $working_dir && $gene_annotation_dir && $gene_annotation_set && $gene_name_type){
  print GREEN, "$usage", RESET;
  exit(1);
}

#Check input dirs
$working_dir = &checkDir('-dir'=>$working_dir, '-clear'=>"no");
$gene_annotation_dir = &checkDir('-dir'=>$gene_annotation_dir, '-clear'=>"no");
$bedtools_bin_dir = &checkDir('-dir'=>$bedtools_bin_dir, '-clear'=>"no");

#Set reference annotation files paths
my $ref_ec_file = "$gene_annotation_dir"."$gene_annotation_set".".ExonContent.bed";
unless (-e $ref_ec_file){
  print RED, "\n\nCould not find ref_ec_file: $ref_ec_file\n\n", RESET;
  exit(1);
}
my $ref_junction_file;
my $outfile;
if ($gene_name_type =~ /symbol/i){
  $ref_junction_file = "$gene_annotation_dir"."$gene_annotation_set".".Genes.junc";
  $outfile = $working_dir . "observed.junctions.anno."."$gene_annotation_set".".tsv";
}elsif($gene_name_type =~ /id/i){
  $ref_junction_file = "$gene_annotation_dir"."$gene_annotation_set".".Genes.gid.junc";
  $outfile = $working_dir . "observed.junctions.anno."."$gene_annotation_set".".gid.tsv";
}else{
  print RED, "\n\nGene name type not recognized: $gene_name_type (must be 'symbol' or 'id')", RESET;
  exit(1);
}
unless (-e $ref_junction_file){
  print RED, "\n\nCould not find ref_junction_file: $ref_junction_file\n\n", RESET;
  exit(1);
}

#Make sure the observed junction file can be found
unless (-e $obs_junction_file){
  print RED, "\n\nObserved junction file ($obs_junction_file) could not be found\n\n", RESET;
  exit(1);
}

#Import the reference junctions (from Ensembl + UCSC + MGC + Refseq + Vega + CCDS)
my %known_junctions;
my %known_donors;
my %known_acceptors;
&importRefJunctions('-ref_junction_file'=>$ref_junction_file);

#Import the reference exons (actually merged exon content blocks (from Ensembl + UCSC + MGC + Refseq + Vega + CCDS)
my %known_ec_blocks;
&importEcBlocks('-ref_ec_file'=>$ref_ec_file);

#Import all observed junctions
my %observed_junctions;
&importObsJunctions('-obs_junction_file'=>$obs_junction_file);

#Determine the exon and splice site skipping of each junction observed (if it known, or anchored to a known splice site at one or both ends)
#&annotateSkipping();
&annotateSkippingBT('-working_dir'=>$working_dir);

#Calculate the junction read count per million junction reads mapped (JPM)
my $grand_count = 0;
foreach my $jid (sort {$observed_junctions{$a}{order} <=> $observed_junctions{$b}{order}} keys %observed_junctions){
  $grand_count += $observed_junctions{$jid}{read_count};
}
foreach my $jid (sort {$observed_junctions{$a}{order} <=> $observed_junctions{$b}{order}} keys %observed_junctions){
  my $jpm = $observed_junctions{$jid}{read_count} * (1000000 / $grand_count);
  $observed_junctions{$jid}{jpm} = $jpm;
}


open(OUT, ">$outfile") || die "\n\nCould not open output file: $outfile\n\n";
print OUT "JID\tRead_Count\tJPM\tIntron_Size\tSplice_Site\tAnchored\tExons_Skipped\tDonors_Skipped\tAcceptors_Skipped\tGene_Name\n";
foreach my $jid (sort {$observed_junctions{$a}{order} <=> $observed_junctions{$b}{order}} keys %observed_junctions){
  my $string = "$jid\t$observed_junctions{$jid}{read_count}\t$observed_junctions{$jid}{jpm}\t$observed_junctions{$jid}{intron_size}\t$observed_junctions{$jid}{splice_site}\t$observed_junctions{$jid}{anchored}\t$observed_junctions{$jid}{exons_skipped}\t$observed_junctions{$jid}{donors_skipped}\t$observed_junctions{$jid}{acceptors_skipped}\t$observed_junctions{$jid}{gid_list}\n";

  #print "$string";
  print OUT "$string";
}
close(OUT);

if ($verbose){ print BLUE, "\n\nPrinted resulting annotated junctions to: $outfile", RESET; }

if ($verbose){ print "\n\n"; }

exit();


######################################################################################################################################
#Import the reference junctions (e.g. Ensembl + UCSC + MGC + Refseq + Vega + CCDS)                                                   #
######################################################################################################################################
sub importRefJunctions{
  my %args = @_;
  my $infile = $args{'-ref_junction_file'};

  if ($verbose){ print BLUE, "\n\nImporting reference junctions", RESET; }
  #Inport the ref junction file.  In this file, the coordinates are always ordered regardless of strand
  #i.e. the 'left' coordinate is always a smaller number than the 'right' coordinate
  #If the strand is '+' then Donor=left and Acceptor=right
  #If the strand is '-' then Donor=right and Acceptor=left
  open(JUNC, "$infile") || die "\n\nCould not open input file: $infile\n\n";
  my $c = 0;
  while(<JUNC>){
    chomp($_);
    my @line = split("\t", $_);
    my $jid = $line[0];
    my $chr = $line[1];
    my $left = $line[2];
    my $right = $line[3];
    my $strand = $line[4];
    my $count = $line[5];
    my $gid_list = $line[6];

    $known_junctions{$jid}=$gid_list;
    $c++;

    if ($strand eq "+"){
      $known_donors{$chr}{$strand}{$left}=$gid_list;
      $known_acceptors{$chr}{$strand}{$right}=$gid_list;

    }elsif($strand eq "-"){
      $known_donors{$chr}{$strand}{$right}=$gid_list;
      $known_acceptors{$chr}{$strand}{$left}=$gid_list;
    }else{
      print RED, "\n\nUnknown strand in reference junctions file", RESET;
      exit(1);
    }
  }
  close(JUNC);

  if ($verbose){ print BLUE, "\n\tImported $c known junctions and associated acceptor/donor pairs\n", RESET; }

  return();
}


######################################################################################################################################
#Import the reference exons (actually merged exon content blocks (from Ensembl + UCSC + MGC + Refseq + Vega + CCDS)                  #
######################################################################################################################################
sub importEcBlocks{
  my %args = @_;
  my $infile = $args{'-ref_ec_file'};

  if ($verbose){ print BLUE, "\n\nImporting reference exon content blocks", RESET; }

  #Import the ref exon content block file.  In this file, the coordinates are always ordered regardless of strand
  #i.e. the 'left' coordinate is always a smaller number than the 'right' coordinate
  open(EXON, "$infile") || die "\n\nCould not open input file: $infile\n\n";
  my $c = 0;
  while(<EXON>){
    chomp($_);
    my @line = split("\t", $_);


    my $chr = $line[0];
    my $left = $line[1];
    my $right = $line[2];
    my $count = $line[3];
    my $strand = $line[4];

    $c++;

    if ($known_ec_blocks{$chr}{$strand}){
      my $ec_ref = $known_ec_blocks{$chr}{$strand};
      $ec_ref->{$c}->{left} = $left;
      $ec_ref->{$c}->{right} = $right;
    }else{
      my %tmp;
      $tmp{$c}{left} = $left;
      $tmp{$c}{right} = $right;
      $known_ec_blocks{$chr}{$strand} = \%tmp;
    }
  }
  close(EXON);
  if ($verbose) { print BLUE, "\n\tImported $c known exon content blocks\n", RESET; }

  return();
}


######################################################################################################################################
#Import the observed junctions from the specified file                                                                               #
######################################################################################################################################
sub importObsJunctions{
  my %args = @_;
  my $infile = $args{'-obs_junction_file'};

  if ($verbose) { print BLUE, "\n\nImporting observed junctions", RESET; }
  #Inport the obs junction file.  In this file, the coordinates are always ordered regardless of strand
  #i.e. the 'left' coordinate is always a smaller number than the 'right' coordinate
  #The strand is unknown...
  open(JUNC, "$infile") || die "\n\nCould not open input file: $infile\n\n";
  my $c = 0;
  my $header = 1;
  my %columns;

  while(<JUNC>){
    chomp($_);
    my @line = split("\t", $_);

    if ($header){
      my $pos = 0;
      foreach my $head (@line){
        $columns{$head}{pos} = $pos;
        $pos++;
      }
      $header = 0;
      next();
    }

    my $jid = $line[$columns{'chr:start-end'}{pos}];
    my $read_count;
    if ($columns{'read_count'}){
      $read_count = $line[$columns{'read_count'}{pos}];
    }elsif($columns{'transcript_count'}){
      $read_count = $line[$columns{'transcript_count'}{pos}];
    }else{
      print RED, "\n\nCould not identify count column\n\n", RESET;
      exit();
    }
    my $intron_size = $line[$columns{'intron_size'}{pos}];
    my $splice_site = $line[$columns{'splice_site'}{pos}];
    
    my $chr;
    my $left;
    my $right;
    my $strand;
    if ($jid =~ /(.*)\:(\d+)\-(\d+)\((.*)\)/){
      $chr = $1;
      $left = $2;
      $right = $3;
      $strand = $4;
    }else{
      print RED, "\n\nObserved junction not understood\n\n", RESET;
      exit(1);
    }

    $observed_junctions{$jid}{read_count} = $read_count;
    $observed_junctions{$jid}{intron_size} = $intron_size;
    $observed_junctions{$jid}{splice_site} = $splice_site;
    $observed_junctions{$jid}{anchored} = "N";
    $observed_junctions{$jid}{gid_list} = "na";
    $c++;

    #First check for an exact match to a known junction (i.e. anchored by both Donor and Acceptor)
    $observed_junctions{$jid}{order} = $c;
    if ($known_junctions{$jid}){
      $observed_junctions{$jid}{anchored} = "DA";
      $observed_junctions{$jid}{gid_list} = $known_junctions{$jid};
    }else{
      #Now check for anchoring on one side only
      #    - If the Donor matches to a gene on the positive strand, the non matching Acceptor should have a higher chromosome coordinate
      #    - If the Donor matches to a gene on the negative strand, the non matching Acceptor should have a lower chromosome coordinate
      #    - If the Acceptor matches to a gene on the positive strand, the non-matching Donor should have a lower chromosome coordinate
      #    - If the Acceptor matches to a gene on the negative strand, the non-matching Donor should have a higher chromosome coordinate
      if ($strand eq "+"){
        if($known_donors{$chr}{'+'}{$left} && $known_acceptors{$chr}{'+'}{$right}){
          $observed_junctions{$jid}{anchored} = "NDA";
          $observed_junctions{$jid}{gid_list} = $known_donors{$chr}{'+'}{$left};
        }elsif($known_donors{$chr}{'+'}{$left}){
          $observed_junctions{$jid}{anchored} = "D";
          $observed_junctions{$jid}{gid_list} = $known_donors{$chr}{'+'}{$left};
        }elsif($known_acceptors{$chr}{'+'}{$right}){
          $observed_junctions{$jid}{anchored} = "A";
          $observed_junctions{$jid}{gid_list} = $known_acceptors{$chr}{'+'}{$right};
        }
      }elsif ($strand eq "-"){
        if($known_donors{$chr}{'-'}{$right} && $known_acceptors{$chr}{'-'}{$left}){
          $observed_junctions{$jid}{anchored} = "NDA";
          $observed_junctions{$jid}{gid_list} = $known_donors{$chr}{'-'}{$right};          
        }elsif ($known_donors{$chr}{'-'}{$right}){
          $observed_junctions{$jid}{anchored} = "D";
          $observed_junctions{$jid}{gid_list} = $known_donors{$chr}{'-'}{$right};
        }elsif($known_acceptors{$chr}{'-'}{$left}){
          $observed_junctions{$jid}{anchored} = "A";
          $observed_junctions{$jid}{gid_list} = $known_acceptors{$chr}{'-'}{$left};
        }
      }
    }

    #Initialize the skipping values
    my $anchored = $observed_junctions{$jid}{anchored};

    if ($anchored eq "N"){
      $observed_junctions{$jid}{exons_skipped} = "na";
      $observed_junctions{$jid}{donors_skipped} = "na";
      $observed_junctions{$jid}{acceptors_skipped} = "na";
    }else{
      $observed_junctions{$jid}{exons_skipped} = 0;
      $observed_junctions{$jid}{donors_skipped} = 0;
      $observed_junctions{$jid}{acceptors_skipped} = 0;
    }

  }
  close(JUNC);

  if ($verbose) { print BLUE, "\n\tImported $c observed junctions\n", RESET; }
  return();
}


##############################################################################################################################################
#Determine the exon and splice site skipping of each junction observed (if it known, or anchored to a known splice site at one or both ends) #
##############################################################################################################################################
sub annotateSkippingBT{
  my %args = @_;
  my $working_dir = $args{'-working_dir'};

  if ($verbose) { print BLUE, "\n\nAnnotating skipping of each junction - exons, then acceptors, then donors - Using BEDTools\n", RESET; }

  #Use BEDTools to determine the exons, acceptors, or donors contained within each exon-exon observed junction (i.e. intron)
  #Do this by writing a temp BED file for each pair of coordinates (of the form: chr, start, end, strand) and then parsing the results file

  #'intersectBed -a exons.bed -b junctions.bed -f 1.0 -s -wa -wb'  
  #The '-f 1.0' option should give the exons that are entirely overlapped by junctions
  #The '-s' option should make this overlap between things on the same strand
  #The '-wa -wb' options, write the original coordinates for exons and junctions (as opposed to the merged coordinates).  Each overlaping pair will be reported as a seperate line
  #Then process the output file and count the exon entries associated with each junction entry
  
  #Example BEDTools command
  #Run BEDTools as follows and 'cut' the column containing the junction ID value.  
  #The occurence count for each of these IDs, should be the number of exons contained within the corresponding junction (i.e. skipped)
  #Unix sort and uniq can even do this counting for us...
  #Then just parse the counts out

  my $temp_obs_junctions = "$working_dir"."ObsJunc.tmp.bed";
  my $temp_known_exons = "$working_dir"."KnownExonContent.tmp.bed";
  my $temp_known_donors = "$working_dir"."KnownDonors.tmp.bed";
  my $temp_known_acceptors = "$working_dir"."KnownAcceptors.tmp.bed";
  my $temp_result = "$working_dir"."Result.tmp.txt";

  #OBSERVED JUNCTIONS
  if ($verbose) { print BLUE, "\n\tLooking for entire exons skipped", RESET; }
  open (TMP_OJ, ">$temp_obs_junctions") || die "\n\nCould not open output file: $temp_obs_junctions";
  #Print print out the observed hmmSplicer junctions (i.e. introns as a temp bed file
  foreach my $jid (sort {$observed_junctions{$a}{order} <=> $observed_junctions{$b}{order}} keys %observed_junctions){
    my $chr;
    my $left;
    my $right;
    my $strand;
    if ($jid =~ /(.*)\:(\d+)\-(\d+)\((.*)\)/){
      $chr = $1;
      $left = $2;
      $right = $3;
      $strand = $4;
    }else{
      print RED, "\n\nObserved junction not understood\n\n", RESET;
      exit(1);
    }
    print TMP_OJ "$chr\t$left\t$right\t$jid\t.\t$strand\n";
  }
  close (TMP_OJ);

  #ENTIRE EXONS
  #Print out all the known exon content blocks
  open (TMP_KE, ">$temp_known_exons") || die "\n\nCould not open output file: $temp_known_exons";
  foreach my $chr (sort keys %known_ec_blocks){
    foreach my $strand (sort keys %{$known_ec_blocks{$chr}}){
      my $ec_ref = $known_ec_blocks{$chr}{$strand};

      foreach my $ec (sort {$ec_ref->{$a}->{left} <=> $ec_ref->{$b}->{left}} keys %{$ec_ref}){
        my $ec_left = $ec_ref->{$ec}->{left};
        my $ec_right = $ec_ref->{$ec}->{right};
        print TMP_KE "$chr\t$ec_left\t$ec_right\tEC\t.\t$strand\n";

      }
    }
  }
  close (TMP_KE);

  my $bed_cmd1 = "$bedtools_bin_dir"."intersectBed -a $temp_known_exons -b $temp_obs_junctions -f 1.0 -s -wa -wb | cut -f 10 | sort | uniq -c > $temp_result";
  if ($verbose) { print BLUE, "\n\t$bed_cmd1", RESET; }
  Genome::Sys->shellcmd(cmd => $bed_cmd1);
  open (COUNTS, "$temp_result") || die "\n\nCould not open temp results file: $temp_result\n\n";
  while(<COUNTS>){
    chomp($_);
    if ($_ =~ /(\d+)\s+(.*)/){
      my $count = $1;
      my $jid = $2;
      if ($observed_junctions{$jid}){
        my $anchored = $observed_junctions{$jid}{anchored};
        unless ($anchored eq "N"){
          $observed_junctions{$jid}{exons_skipped} = $count;
        }
      }else{
        print RED, "\n\nUnrecognized junction id: $jid\n\n", RESET;
        exit(1);
      }
    }else{
      print RED, "\n\nEntry in results file not understood: $_\n\n", RESET;
      exit(1);
    }
  }
  close(COUNTS);


  #DONORS
  if ($verbose){ print BLUE, "\n\tLooking for donors skipped", RESET; }
  open (TMP_KD, ">$temp_known_donors") || die "\n\nCould not open output file: $temp_known_donors";
  foreach my $chr (sort keys %known_donors){
    foreach my $strand (sort keys %{$known_donors{$chr}}){
      foreach my $donor (sort keys %{$known_donors{$chr}{$strand}}){
        my $donor_p = $donor+1;
        print TMP_KD "$chr\t$donor\t$donor_p\tD\t.\t$strand\n";
      }
    }
  }
  close (TMP_KD);

  my $bed_cmd2 = "$bedtools_bin_dir"."intersectBed -a $temp_known_donors -b $temp_obs_junctions -f 1.0 -s -wa -wb | cut -f 10 | sort | uniq -c > $temp_result";
  if ($verbose){ print BLUE, "\n\t$bed_cmd2", RESET; }
  Genome::Sys->shellcmd(cmd => $bed_cmd2);
  open (COUNTS, "$temp_result") || die "\n\nCould not open temp results file: $temp_result\n\n";
  while(<COUNTS>){
    chomp($_);
    if ($_ =~ /(\d+)\s+(.*)/){
      my $count = $1;
      my $jid = $2;
      if ($observed_junctions{$jid}){
        my $anchored = $observed_junctions{$jid}{anchored};
        unless ($anchored eq "N"){
          $observed_junctions{$jid}{donors_skipped} = $count;
        }
      }else{
        print RED, "\n\nUnrecognized junction id: $jid\n\n", RESET;
        exit(1);
      }
    }else{
      print RED, "\n\nEntry in results file not understood: $_\n\n", RESET;
      exit(1);
    }
  }
  close(COUNTS);


  #ACCEPTORS
  if ($verbose){ print BLUE, "\n\tLooking for acceptors skipped", RESET; }
  open (TMP_KA, ">$temp_known_acceptors") || die "\n\nCould not open output file: $temp_known_acceptors";
  foreach my $chr (sort keys %known_acceptors){
    foreach my $strand (sort keys %{$known_acceptors{$chr}}){
      foreach my $acceptor (sort keys %{$known_acceptors{$chr}{$strand}}){
        my $acceptor_p = $acceptor+1;
        print TMP_KA "$chr\t$acceptor\t$acceptor_p\tA\t.\t$strand\n";
      }
    }
  }
  close (TMP_KA);

  my $bed_cmd3 = "$bedtools_bin_dir"."intersectBed -a $temp_known_acceptors -b $temp_obs_junctions -f 1.0 -s -wa -wb | cut -f 10 | sort | uniq -c > $temp_result";
  if ($verbose){ print BLUE, "\n\t$bed_cmd3", RESET; }
  Genome::Sys->shellcmd(cmd => $bed_cmd3);
  open (COUNTS, "$temp_result") || die "\n\nCould not open temp results file: $temp_result\n\n";
  while(<COUNTS>){
    chomp($_);
    if ($_ =~ /(\d+)\s+(.*)/){
      my $count = $1;
      my $jid = $2;
      if ($observed_junctions{$jid}){
        my $anchored = $observed_junctions{$jid}{anchored};
        unless ($anchored eq "N"){
          $observed_junctions{$jid}{acceptors_skipped} = $count;
        }
      }else{
        print RED, "\n\nUnrecognized junction id: $jid\n\n", RESET;
        exit(1);
      }
    }else{
      print RED, "\n\nEntry in results file not understood: $_\n\n", RESET;
      exit(1);
    }
  }
  close(COUNTS);



  
  #Clean up the temp files
  my $rm_cmd = "rm -f $temp_obs_junctions $temp_known_exons $temp_known_donors $temp_known_acceptors $temp_result";
  if ($verbose){ print BLUE, "\n\nCleaning up ...\n$rm_cmd", RESET; }
  Genome::Sys->shellcmd(cmd => $rm_cmd);

  return();
}


