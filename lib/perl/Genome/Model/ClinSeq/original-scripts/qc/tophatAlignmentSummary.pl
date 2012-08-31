#!/usr/bin/env genome-perl
#Written by Malachi Griffith

#This script will examine a tophat alignment results directory and generate basic statistics on the alignmnents
#As a tool it should take in input files and produce output files

#Main inputs:
#RNA-seq Build Dir (e.g. /gscmnt/gc8001/info/model_data/2881643231/build117377906/)
#Annotation Build Dir (e.g. /gscmnt/ams1102/info/model_data/2771411739/build106409619/)

#Load modules
use strict;
use warnings;
use Getopt::Long;
use Term::ANSIColor qw(:constants);
use Data::Dumper;
use above 'Genome'; #remove the 'above' when we turn this into a module
use Genome::Model::ClinSeq::Util qw(:all);

my $script_dir;
use Cwd 'abs_path';
BEGIN{
  if (abs_path($0) =~ /(.*\/).*\/.*\.pl/){
    $script_dir = $1;
  }
}

#Input parameters
my $reference_fasta_file = '';
my $tophat_alignment_dir = '';
my $reference_annotations_dir = '';
my $working_dir = '';
my $verbose = 0;
my $clean = 0;

GetOptions ('reference_fasta_file=s'=>\$reference_fasta_file, 'tophat_alignment_dir=s'=>\$tophat_alignment_dir, 'reference_annotations_dir=s'=>\$reference_annotations_dir,
 	    'working_dir=s'=>\$working_dir, 'verbose=i'=>\$verbose, 'clean=i'=>\$clean);


my $usage=<<INFO;

  Example usage: 
  
  tophatAlignmentStats.pl  --reference_fasta_file='/gscmnt/sata420/info/model_data/2857786885/build102671028/all_sequences.fa'  --tophat_alignment_dir='/gscmnt/gc2014/info/model_data/2880794541/build115909743/alignments/'  --reference_annotations_dir='/gscmnt/sata132/techd/mgriffit/reference_annotations/hg19/'  --working_dir='/gscmnt/sata132/techd/mgriffit/hgs/hg1/qc/tophat/'
  
  Intro:
  This script summarizes results from a tophat alignment directory and writes resulting stats and figures to a working directory

  Details:
  --reference_fasta_file          Reference fasta file that was used for Tophat mapping
  --tophat_alignment_dir          The 'alignment' dir created by a Tophat run
  --reference_annotations_dir     Directory containing the reference junctions to be compared against
                                  For example: /gscmnt/sata132/techd/mgriffit/reference_annotations/hg19/ALL.Genes.junc
  --working_dir                   Directory where results will be stored
  --verbose                       To display more output, set to 1
  --clean                         To clobber the top dir and create everything from scratch, set to 1

INFO

unless ($reference_fasta_file && $tophat_alignment_dir && $reference_annotations_dir && $working_dir){
  print GREEN, "$usage", RESET;
  exit(1);
}

#Check input directories and files
$tophat_alignment_dir = &checkDir('-dir'=>$tophat_alignment_dir, '-clear'=>"no");
if ($clean){
  $working_dir = &checkDir('-dir'=>$working_dir, '-clear'=>"yes", '-recursive'=>"yes");
}else{
  $working_dir = &checkDir('-dir'=>$working_dir, '-clear'=>"no");
}
unless (-e $reference_annotations_dir){
  print RED, "\n\nCould not find reference junctions file: $reference_annotations_dir\n\n", RESET;
  exit(1);
}
#Since UCSC annotations are used for known splice site information, a separate reference fasta file will need to be specified
#The other one used for the alignments and specified in the command can come from the TGI API
my $annotation_reference_fasta_file = "$reference_annotations_dir/ref_fasta/all_sequences.fa";
unless (-e $annotation_reference_fasta_file){
  print RED, "\n\nCould not find annotation reference fasta: $annotation_reference_fasta_file\n\n", RESET;
  exit();
}

my $tophat_stats_file = $tophat_alignment_dir . "alignment_stats.txt";
my $tophat_junctions_bed_file = $tophat_alignment_dir . "junctions.bed";
my $new_tophat_junctions_bed_file = $working_dir . "observed.junctions.bed";
my $tophat_junctions_junc_file = $working_dir . "observed.junctions.junc";
my $tophat_junctions_anno_file = $working_dir . "observed.junctions.splicesites.junc";
my @gene_name_type = qw ( symbol id );
my @gene_annotation_sets = qw (ALL Ensembl);

#Make a copy of the junctions file and alignment stats file
my $cp_cmd1 = "cp $tophat_junctions_bed_file $new_tophat_junctions_bed_file";
if ($verbose){ print YELLOW, "\n\n$cp_cmd1", RESET; }
Genome::Sys->shellcmd(cmd => $cp_cmd1);

my $cp_cmd2 = "cp $tophat_stats_file $working_dir";
if ($verbose){ print YELLOW, "\n\n$cp_cmd2", RESET; }
Genome::Sys->shellcmd(cmd => $cp_cmd2);


#Convert observed junctions.bed to an observed .junc file
my $bed_to_junc_cmd = "cat $new_tophat_junctions_bed_file | "."$script_dir"."misc/bed2junc.pl > $tophat_junctions_junc_file";
if ($verbose){ print YELLOW, "\n\n$bed_to_junc_cmd", RESET; }
Genome::Sys->shellcmd(cmd => $bed_to_junc_cmd);


#Go through the .junc file and infer the splice site of each observed junction
&inferSpliceSite('-infile'=>$tophat_junctions_junc_file, '-outfile'=>$tophat_junctions_anno_file, '-reference_fasta_file'=>$reference_fasta_file);


#Grab the known junctions file for each annotation set and infer the splice sites of these known junctions
foreach my $gene_annotation_set (@gene_annotation_sets){
  my $ref_junction_file = "$reference_annotations_dir"."$gene_annotation_set".".Genes.junc"; 
  my $known_junctions_file = $working_dir . "$gene_annotation_set".".junctions.junc";
  my $known_junctions_anno_file = $working_dir . "$gene_annotation_set".".junctions.splicesites.junc";
  #Make a simple junctions file from the source annotation file
  #chr:start-end	transcript_count (where 'transcript_count' is the number of known transcripts in the file that use that junction)
  open (IN, "$ref_junction_file") || die "\n\nCould not open input file: $ref_junction_file\n\n";
  open (OUT, ">$known_junctions_file") || die "\n\nCould not open output file: $known_junctions_file\n\n";
  print OUT "chr:start-end\ttranscript_count\n";
  while(<IN>){
    chomp($_);
    my @line = split("\t", $_);
    print OUT "$line[0]\t$line[5]\n";
  }
  close(IN);
  close(OUT);
  &inferSpliceSite('-infile'=>$known_junctions_file, '-outfile'=>$known_junctions_anno_file, '-reference_fasta_file'=>$annotation_reference_fasta_file);
}


#Now annotate observed exon-exon junctions against databases of known junctions
foreach my $gene_name_type (@gene_name_type){
  foreach my $gene_annotation_set (@gene_annotation_sets){
    my $cmd = "$script_dir/rnaseq/annotateObservedJunctions.pl  --obs_junction_file=$tophat_junctions_anno_file  --bedtools_bin_dir='/gsc/bin/'  --working_dir=$working_dir  --gene_annotation_dir=$reference_annotations_dir  --gene_annotation_set='$gene_annotation_set'  --gene_name_type='$gene_name_type'  --verbose=$verbose";
    if ($verbose){ print YELLOW, "\n\n$cmd", RESET };
    Genome::Sys->shellcmd(cmd => $cmd);
  }
}


#Ensembl based gene/transcript level expression calculations:
my $ensembl_observed_junction_file = $working_dir . "observed.junctions.anno.Ensembl.gid.tsv";
my $ensembl_junction_file = "$reference_annotations_dir"."Ensembl.Genes.junc"; 
my $ensembl_junction_gid_file = "$reference_annotations_dir"."Ensembl.Genes.gid.junc";
my $ensembl_gene_info_file = "$reference_annotations_dir"."transcript_to_gene/Ensembl.Genes.info";
&calculateGeneExpression('-working_dir'=>$working_dir, '-known_junction_file'=>$ensembl_junction_file, '-known_junction_gid_file'=>$ensembl_junction_gid_file, '-observed_junction_file'=>$ensembl_observed_junction_file, '-gene_info_file'=>$ensembl_gene_info_file);


#Feed resulting files into R and generate statistics (also supply known junctions file used for the analysis)
#- Summarize the alignment stats file using an R script
#- Basic junction stats: 
#  - Total junctions observed 
#  - Total known junctions observed
#  - Proportion of all known junctions observed
#  - Total exon skipping junctions observed (and proportion of the library)
#  - Total novel exon skipping junctions observed (and proportion of the library)
#- Pie chart of splice sites observed (GC-AG, GC-AG, etc.) - compare to a pie chart of all known splice sites ...
#- Pie chart of anchor types (DA, NDA, D, A, N) - Number and Percentage of reads corresponding to each type
#- Percentage of all junction mapping reads corresponding to known junctions
#- Distribution of exon-exon junction read counts
#- Display read count distribution at both the gene and junction level
#- Expression distribution bias.  Percentage of all reads consumed by top N .. M % of junctions/genes
#- For exon-skipping events, display the proportion that are 1S, 2S, 3S, etc. - Show for all known, known observed, novel observed
#- How many genes are covered over the majority of their junctions (25%, 50%, 75%, 90%, 95%, 100%)?
#- Produce ranked gene expression lists based on exon-junction values
#- What is the distribution of observed intron sizes versus known intron sizes
#- Produce a Top N% expressed file

my $r_cmd = "$script_dir/qc/tophatAlignmentSummary.R $working_dir";
my $r_cmd_stdout = "$working_dir"."tophatAlignmentSummary.R.stdout";
my $r_cmd_stderr = "$working_dir"."tophatAlignmentSummary.R.stderr";
if ($verbose){ 
  print YELLOW, "\n\nRunning: $r_cmd", RESET;
}else{
  $r_cmd .= " 1>$r_cmd_stdout 2>$r_cmd_stderr";
}
Genome::Sys->shellcmd(cmd => $r_cmd);


#Get some basic info from the Tophat stats file and append to the existing Stats.tsv file
#- Total and percent of total for the following: reads, reads mapped, unmapped reads, multiple hit reads
open (TOPHAT, "$tophat_stats_file") || die "\n\nCould not open tophat stats file\n\n";
my $total_reads=1;
my $unmapped_reads=0;
my $multimap_reads=0;
my $mapped_reads=0;
my $percent_mapped=0;
my $percent_unmapped=0;
my $percent_multimap=0;
while(<TOPHAT>){
  chomp($_);
  unless ($_ =~ /^\#/){next();}
  if ($_ =~ /Total\s+Reads\:\s+(\d+)/){$total_reads=$1;}
  if ($_ =~ /Unmapped\s+Reads\:\s+(\d+)/){$unmapped_reads=$1;}
  if ($_ =~ /Multiple\s+Hit\s+Reads\:\s+(\d+)/){$multimap_reads=$1;}
  if ($_ =~ /Total\s+Reads\s+Mapped\:\s+(\d+)/){$mapped_reads=$1;}
  $percent_mapped = sprintf("%.2f", (($mapped_reads/$total_reads)*100));
  $percent_unmapped = sprintf("%.2f", (($unmapped_reads/$total_reads)*100));
  $percent_multimap = sprintf("%.2f", (($multimap_reads/$total_reads)*100));
}
close(TOPHAT);

my $new_stats_file = $working_dir . "summary/Stats.tsv";
open (STATS, ">>$new_stats_file") || die "\n\nCould not open new Stats.tsv file\n\n";
print STATS "Total Reads\t$total_reads\tRNA-seq\tAlignments\tCount\tTotal reads used for Tophat alignments\n";
print STATS "Mapped Reads\t$mapped_reads\tRNA-seq\tAlignments\tCount\tTotal reads with at least one mapping by Tophat\n";
print STATS "Percent Mapped Reads\t$percent_mapped\tRNA-seq\tAlignments\tPercent\tPercent of reads with at least one mapping by Tophat\n";
print STATS "MultiMap Reads\t$multimap_reads\tRNA-seq\tAlignments\tCount\tTotal reads with more than one mapping by Tophat\n";
print STATS "Percent MultiMap Reads\t$percent_multimap\tRNA-seq\tAlignments\tPercent\tPercent of reads with more than one mapping by Tophat\n";
print STATS "Unmapped Reads\t$unmapped_reads\tRNA-seq\tAlignments\tCount\tTotal reads with NO mappings by Tophat\n";
print STATS "Percent Unmapped Reads\t$percent_unmapped\tRNA-seq\tAlignments\tPercent\tPercent of reads with NO mappings by Tophat\n";
close(STATS);

if ($verbose){print "\n\n"};

exit();



############################################################################################################################
#Infer strand from alignment to reference genome                                                                           #
############################################################################################################################
sub inferSpliceSite{
  my %args = @_; 
  my $infile = $args{'-infile'};
  my $outfile = $args{'-outfile'};
  my $reference_fasta_file = $args{'-reference_fasta_file'};

  #Load in the junctions from the input file
  my %junctions;
  open(JUNC, "$infile") || die "\n\nCould not open input file: $infile\n\n";
  my $header = 1;
  my $header_line;
  my %columns;
  my $o = 0;
  while(<JUNC>){
    chomp($_);
    my @line = split("\t", $_);
    if ($header){
      $header_line = $_;
      my $p = 0;
      foreach my $col (@line){
        $columns{$col}{position} = $p;
        $p++;
      }
      $header = 0;
      next();
    }
    $o++;
    my $jid = $line[$columns{'chr:start-end'}{position}];
    my $read_count;
    if ($columns{'read_count'}){
      $read_count = $line[$columns{'read_count'}{position}];
    }elsif($columns{'transcript_count'}){
      $read_count = $line[$columns{'transcript_count'}{position}];
    }else{
      print RED, "\n\nCould not identify count column\n\n", RESET;
      exit();
    }
    $junctions{$jid}{read_count} = $read_count;
    $junctions{$jid}{order} = $o;
  }
  close(JUNC);

  #Determine the splice site by comparison back to the reference genome...
  if ($verbose){print BLUE, "\n\nAttempting to determine splice site for each junction in:\n\t$infile", RESET;}

  open (REF, "$reference_fasta_file") || die "\n\nCould not open fasta file: $reference_fasta_file";
  my $tmp = $/;
  $/ = "\n>";  # read by FASTA record

  while (<REF>){
    chomp $_;
    my $chr_seq = $_;
    my ($chr) = $chr_seq =~ /^>*(\S+)/;  # parse ID as first word in FASTA header

    #Add 'chr' to the begining of the chr name only if it isnt there already
    unless ($chr =~ /^chr/){
      $chr = "chr".$chr;
    }
    $chr_seq =~ s/^>*.+\n//;  # remove FASTA header
    $chr_seq =~ s/\n//g;  # remove endlines

    my $chr_length = length($_);
    if ($verbose){print BLUE, "\n\tFound $chr sequence (length = $chr_length)", RESET;}

    #Now go through the junctions found for this chromosome and look for donor/acceptor splice sites at the coordinates reported
    #SPLICE_SITES = ["GT-AG", "CT-AC", "GC-AG", "CT-GC", "AT-AC", "GT-AT"]

    foreach my $j (sort keys %junctions){
      if ($j =~ /(.*)\:(\d+)\-(\d+)/){
        my $j_chr = $1;
        my $left = $2;
        my $right = $3;
        unless($chr eq $j_chr){
          next();
        }
        my $intron_size = ($right - $left)+1;
        $junctions{$j}{intron_size} = $intron_size;
        my $left_dn = uc(substr($chr_seq, $left, 2));
        my $right_dn = uc(substr($chr_seq, $right-3, 2));

        #print "\n\t\tDEBUG: $left_dn ... $right_dn";
        #Strand is assigned by Tophat...
        if ($left_dn eq "GT" && $right_dn eq "AG"){
          $junctions{$j}{splice_site} = "GT-AG";
        }elsif($left_dn eq "CT" && $right_dn eq "AC"){
          $junctions{$j}{splice_site} = "GT-AG";
        }elsif($left_dn eq "GC" && $right_dn eq "AG"){
          $junctions{$j}{splice_site} = "GC-AG";
        }elsif($left_dn eq "CT" && $right_dn eq "GC"){
          $junctions{$j}{splice_site} = "GC-AG";
        }elsif($left_dn eq "AT" && $right_dn eq "AC"){
          $junctions{$j}{splice_site} = "AT-AC";     
        }elsif($left_dn eq "GT" && $right_dn eq "AT"){
          $junctions{$j}{splice_site} = "AT-AC";
        }else{
          $junctions{$j}{splice_site} = "Other";
        }
      }else{
        print RED, "\n\nObserved junction not understood\n\n", RESET;
        exit();
      }
      #print YELLOW, "\n\t$junctions{$j}{strand}\t$junctions{$j}{splice_site}", RESET;
    }
  }
  close(REF);
  $/ = $tmp;

  #Print out the strand inferred junctions to a new file
  open (OUT, ">$outfile") || die "\n\nCould not open output file: $outfile\n\n";
  print OUT "$header_line\tintron_size\tsplice_site\n";
  foreach my $j (sort {$junctions{$a}{order} <=> $junctions{$b}{order}} keys %junctions){

    unless (defined($junctions{$j}{read_count}) && defined($junctions{$j}{intron_size}) && defined($junctions{$j}{splice_site})){
      print RED, "\n\nUndefined value in output!\n\n", RESET;
      print "$header_line\tintron_size\tsplice_site\n";
      print "$j\t$junctions{$j}{read_count}\t$junctions{$j}{intron_size}\t$junctions{$j}{splice_site}\n";
      exit();
    }
    print OUT "$j\t$junctions{$j}{read_count}\t$junctions{$j}{intron_size}\t$junctions{$j}{splice_site}\n";

  }
  close(OUT);

  return();
}


############################################################################################################################
#Ensembl based gene/transcript level expression calculations:                                                              #
############################################################################################################################
sub calculateGeneExpression{
  my %args = @_;
  my $working_dir = $args{'-working_dir'};
  my $known_junction_file = $args{'-known_junction_file'};
  my $known_junction_gid_file = $args{'-known_junction_gid_file'};
  my $observed_junction_file = $args{'-observed_junction_file'};
  my $gene_info_file = $args{'-gene_info_file'};

  #Check for existance of files
  unless (-e $known_junction_file){
    print RED, "\n\nKnown junction file not found: $known_junction_file\n\n", RESET;
    exit();
  }
  unless (-e $known_junction_gid_file){
    print RED, "\n\nKnown junction gid file not found: $known_junction_gid_file\n\n", RESET;
    exit();
  }
  unless (-e $observed_junction_file){
    print RED, "\n\nObserved junction file not found: $observed_junction_file\n\n", RESET;
    exit();
  }
  unless (-e $gene_info_file){
    print RED, "\n\nGene info file not found: $gene_info_file\n\n", RESET;
    exit();
  }



  my $entrez_ensembl_data = &loadEntrezEnsemblData();


  #Calculate gene/transcript level read counts and expression estimates from exon-exon junction counts
  #Gene/transcript level expression = (sum of exon junction read counts for a gene / number of exon-exon junctions of that gene) => then normalized to per million junction mapped reads (JPJM)
  #Only known 'DA' junctions among the observed junctions will be used for these calculations
  #Store the proportion of junctions of the gene that were observed at least 1X, 5X, 10X, etc.
  #Make sure the output file has both ENSG/ENST ID, Gene Symbol, and Mapped Gene Symbol
  #Also output the number of known exon-exon junctions of the gene
  #Note that one junction can correspond to multiple genes - in this implementation reads for these junctions can be counted multiple times (once for each gene they correspond to)??
  #Create a gene/transcript expression record for every gene/transcript that has at least one exon-exon junction regardless of whether it was detected or not

  #1.) Build a map of:
  #-   all genes to gene names and list of transcripts
  #-   all transcripts to gene names and gene ids
  #-   also store a list of known junctions keyed on jid for convenient lookups
  my %genes;
  my %transcripts;
  my %known_junctions;
  my %observed_known_junctions;
  if ($verbose){print BLUE, "\n\nImporting gene, gene name, and transcript mappings from: $gene_info_file", RESET;}
  open (GENE, "$gene_info_file") || die "\n\nCould not open gene info file: $gene_info_file\n\n";
  while(<GENE>){
    chomp($_);
    if ($_ =~ /^\#/){
      next();
    }
    my @line = split("\t", $_);
    my $enst_id = $line[0];
    my $ensg_id = $line[1];
    my $name = $line[3];
    $transcripts{$enst_id}{ensg_id} = $ensg_id;
    $transcripts{$enst_id}{name} = $name;
    $genes{$ensg_id}{ensg_id} = $ensg_id;
    $genes{$ensg_id}{name} = $name;
  }
  close(GENE);

  #2.) Build a map of all known junctions to their: gene_ids, transcript_ids.  i.e.
  #    - Add junctions to transcript hash create above -> also store known junction count
  if ($verbose){print BLUE, "\n\nImporting known junction to transcript mappings from: $known_junction_file", RESET;}
  open (JUNC1, "$known_junction_file") || die "\n\nCould not open known junction file: $known_junction_file\n\n";
  while(<JUNC1>){
    chomp($_);
    my @line = split("\t", $_);
    my $j = $line[0];
    my $chr = $line[1];
    my $count = $line[5];
    my $tid_string = $line[7];
    my @tids = split(",", $tid_string);

    $known_junctions{$j}=1;

    #Attach junction list to each transcript record
    foreach my $enst_id (@tids){
      if (defined($transcripts{$enst_id}{junction_list})){
        my $jlist = $transcripts{$enst_id}{junction_list};
        $jlist->{$j}->{count} = $count;
      }else{
        my %jlist;
        $jlist{$j}{count} = $count;
        $transcripts{$enst_id}{junction_list} = \%jlist;
        $transcripts{$enst_id}{chromosome} = $chr;
      }
    }
  }
  close(JUNC1);

  #    - Add junctions to gene hash created above -> also store known junction count
  if ($verbose){print BLUE, "\n\nImporting known junction to gene mappings from: $known_junction_gid_file", RESET;}
  open (JUNC2, "$known_junction_gid_file") || die "\n\nCould not open known junction file: $known_junction_gid_file\n\n";
  while(<JUNC2>){
    chomp($_);
    my @line = split("\t", $_);
    my $j = $line[0];
    my $chr = $line[1];
    my $count = $line[5];
    my $gid_string = $line[6];
    my @gids = split(",", $gid_string);

    #Attach junction list to each gene record
    foreach my $ensg_id (@gids){
      if (defined($genes{$ensg_id}{junction_list})){
        my $jlist = $genes{$ensg_id}{junction_list};
        $jlist->{$j}->{count} = $count;
      }else{
        my %jlist;
        $jlist{$j}{count} = $count;
        $genes{$ensg_id}{junction_list} = \%jlist;
        $genes{$ensg_id}{chromosome} = $chr;
      }
    }
  }
  close(JUNC2);

  #3.) Build a hash of all observed junctions (only store 'DA' junctions) and their read counts
  #    - Store the grand read count of the library
  #    - Store the grand number of junctions observed for the library
  my $grand_observed_junction_count = 0;
  my $grand_read_count = 0;
  if ($verbose){print BLUE, "\n\nImporting observed junctions and their counts from: $observed_junction_file", RESET;}
  open (JUNC3, "$observed_junction_file") || die "\n\nCould not open observed junction file: $observed_junction_file\n\n";
  while(<JUNC3>){
    chomp($_);
    if ($_ =~ /^JID/){
      next();
    }
    my @line = split("\t", $_);
    my $jid = $line[0];
    my $read_count = $line[1];

    #Only store know junctions ('DA' anchored)
    unless($known_junctions{$jid}){
      next();
    }
    $observed_known_junctions{$jid}{read_count} = $read_count;
    $grand_read_count += $read_count;
  }
  close(JUNC3);
  $grand_observed_junction_count = keys %observed_known_junctions;
  if ($verbose){print BLUE "\n\tFound $grand_observed_junction_count observed known junctions with a grand read count of $grand_read_count", RESET;}


  #4.) Now go through each gene/transcript in the gene/transcript hash, and go through each junction associated with it
  #    - Count the number of junctions actually observed - remove genes/transcripts with 0 known junctions
  #    - Store the cumulative read count for this gene/transcript
  #    - Store the proportion of junctions of the gene that were observed at least 1X, 5X, 10X, etc.
  #    - Gene/transcript level expression = (sum of exon junction read counts for a gene / number of exon-exon junctions of that gene) => then normalized to per million junction mapped reads (JPJM)
  if ($verbose){print BLUE, "\n\nCalculating gene/transcript level expression values", RESET;}
  my $genes_ref = \%genes;
  my $transcripts_ref = \%transcripts;
  my %lists;
  $lists{'transcripts'}{features} = $transcripts_ref;
  $lists{'transcripts'}{outfile} = $working_dir ."Ensembl.Junction.TranscriptExpression.tsv";
  $lists{'genes'}{features} = $genes_ref;
  $lists{'genes'}{outfile} = $working_dir ."Ensembl.Junction.GeneExpression.tsv";

  foreach my $feature_type (keys %lists){
    my $features = $lists{$feature_type}{features};
    my $feature_count = keys %{$features};
    if ($verbose){print BLUE, "\n\nProcessing $feature_count features - removing features with no junctions", RESET;}

    foreach my $fid (keys %{$features}){
      my $junction_list = $features->{$fid}->{junction_list};
      my $known_junction_count = keys %{$junction_list};
      
      #print "\n\tknown_junction count: $known_junction_count";
      if ($known_junction_count > 0){
        $features->{$fid}->{known_junction_count} = $known_junction_count;

        #values to calculate
        my $feature_read_count = 0;
        my $junctions_1x = 0;
        my $junctions_2x = 0;
        my $junctions_5x = 0;
        my $junctions_10x = 0;
        my $junctions_20x = 0;
        my $junctions_50x = 0;
        my $junctions_100x = 0;
        my $junctions_500x = 0;
        my $junctions_1000x = 0;
        foreach my $jid (keys %{$junction_list}){
          my $read_count = 0;
          if ($observed_known_junctions{$jid}){
            $read_count = $observed_known_junctions{$jid}{read_count};
            $feature_read_count += $read_count;
            if ($read_count >= 1){$junctions_1x++;}
            if ($read_count >= 2){$junctions_2x++;}
            if ($read_count >= 5){$junctions_5x++;}
            if ($read_count >= 10){$junctions_10x++;}
            if ($read_count >= 20){$junctions_20x++;}
            if ($read_count >= 50){$junctions_50x++;}
            if ($read_count >= 100){$junctions_100x++;}
            if ($read_count >= 500){$junctions_500x++;}
            if ($read_count >= 1000){$junctions_1000x++;}
          }
        }
        my $junctions_1x_p = sprintf ("%.2f", (($junctions_1x/$known_junction_count)*100));
        my $junctions_2x_p = sprintf ("%.2f", (($junctions_2x/$known_junction_count)*100));
        my $junctions_5x_p = sprintf ("%.2f", (($junctions_5x/$known_junction_count)*100));
        my $junctions_10x_p = sprintf ("%.2f", (($junctions_10x/$known_junction_count)*100));
        my $junctions_20x_p = sprintf ("%.2f", (($junctions_20x/$known_junction_count)*100));
        my $junctions_50x_p = sprintf ("%.2f", (($junctions_50x/$known_junction_count)*100));
        my $junctions_100x_p = sprintf ("%.2f", (($junctions_100x/$known_junction_count)*100));
        my $junctions_500x_p = sprintf ("%.2f", (($junctions_500x/$known_junction_count)*100));
        my $junctions_1000x_p = sprintf ("%.2f", (($junctions_1000x/$known_junction_count)*100));

        #calculate expression level
        #    - Gene/transcript level expression = (sum of exon junction read counts for a gene / number of exon-exon junctions of that gene) => then normalized to per million junction mapped reads (JPJM)
        my $jpj = ($feature_read_count / $known_junction_count); #Junction count normalized to number of junctions for the feature
        my $jpjm = $jpj * (1000000/$grand_read_count); #Further normalized to be per million junction mapping reads

        #Store values for this feature
        $features->{$fid}->{read_count} = $feature_read_count;
        $features->{$fid}->{jpj} = $jpj;
        $features->{$fid}->{jpjm} = $jpjm;
        $features->{$fid}->{junctions_1x_p} = $junctions_1x_p;
        $features->{$fid}->{junctions_2x_p} = $junctions_2x_p;
        $features->{$fid}->{junctions_5x_p} = $junctions_5x_p;
        $features->{$fid}->{junctions_10x_p} = $junctions_10x_p;
        $features->{$fid}->{junctions_20x_p} = $junctions_20x_p;
        $features->{$fid}->{junctions_50x_p} = $junctions_50x_p;
        $features->{$fid}->{junctions_100x_p} = $junctions_100x_p;
        $features->{$fid}->{junctions_500x_p} = $junctions_500x_p;
        $features->{$fid}->{junctions_1000x_p} = $junctions_1000x_p;
      }else{
        delete $features->{$fid};
      }
    }
  }

  #5.) Create output files, one for gene-level and one for transcript-level
  if ($verbose){print BLUE, "\n\nCreating output files", RESET;}
  foreach my $feature_type (keys %lists){
    my $features = $lists{$feature_type}{features};
    my $outfile = $lists{$feature_type}{outfile};
    my $feature_count = keys %{$features};

    if ($verbose){print BLUE, "\n\nProcessing $feature_count features and printing to: $outfile", RESET;}
    open (OUT, ">$outfile") || die "\n\nCould not open output file: $outfile\n\n", RESET;
    print OUT "fid\tensg_id\tgene_name\tmapped_gene_name\tchromosome\tknown_junction_count\tread_count\tjpjm\tjunctions_1x_p\tjunctions_2x_p\tjunctions_5x_p\tjunctions_10x_p\tjunctions_20x_p\tjunctions_50x_p\tjunctions_100x_p\tjunctions_500x_p\tjunctions_1000x_p\n";
    foreach my $fid (sort {$features->{$b}->{jpjm} <=> $features->{$a}->{jpjm}} keys %{$features}){

      #Attempt to map current gene name to an entrez gene name
      my $gene_name = $features->{$fid}->{name};
      my $mapped_gene_name = &fixGeneName('-gene'=>$gene_name, '-entrez_ensembl_data'=>$entrez_ensembl_data, '-verbose'=>0);

      print OUT "$fid\t$features->{$fid}->{ensg_id}\t$gene_name\t$mapped_gene_name\t$features->{$fid}->{chromosome}\t$features->{$fid}->{known_junction_count}\t$features->{$fid}->{read_count}\t$features->{$fid}->{jpjm}\t$features->{$fid}->{junctions_1x_p}\t$features->{$fid}->{junctions_2x_p}\t$features->{$fid}->{junctions_5x_p}\t$features->{$fid}->{junctions_10x_p}\t$features->{$fid}->{junctions_20x_p}\t$features->{$fid}->{junctions_50x_p}\t$features->{$fid}->{junctions_100x_p}\t$features->{$fid}->{junctions_500x_p}\t$features->{$fid}->{junctions_1000x_p}\n";
    }
    close(OUT);

  }
  #print Dumper %transcripts;
  #print Dumper %genes;

  return();
}


