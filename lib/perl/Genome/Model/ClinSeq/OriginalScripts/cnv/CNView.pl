#!/usr/bin/env genome-perl
#Written by Malachi Griffith

use strict;
use Getopt::Long;
use Term::ANSIColor qw(:constants);
use Data::Dumper;
use File::Basename;
use above 'Genome';
use Genome::Model::ClinSeq::Util qw(:all);

my $script_dir;
use Cwd 'abs_path';
BEGIN{
  if (abs_path($0) =~ /(.*\/)(.*\/).*\.pl/){
    $script_dir = $1.$2;
  }
}

#Required
my $reference_build = '';
my $cnv_file = '';
my $model_build = '';
my $working_dir = '';
my $sample_name = '';

#Optional
my $gene_targets_file = '';
my $name = '';
my $cnv_results_file = '';
my $ideogram_file = '';
my $verbose = '';
my $chr = '';
my $chr_start = '';
my $chr_end = '';
my $image_type = '';
my $force = '';

GetOptions ('reference_build=s'=>\$reference_build, 'cnv_file=s'=>\$cnv_file, 'working_dir=s'=>\$working_dir, 'sample_name=s'=>\$sample_name, 'model_build=s'=>\$model_build,
            'gene_targets_file=s'=>\$gene_targets_file, 'name=s'=>\$name, 'cnv_results_file=s'=>\$cnv_results_file, 'ideogram_file=s'=>\$ideogram_file, 'verbose=s'=>\$verbose,
            'chr=s'=>\$chr, 'chr_start=i'=>\$chr_start, 'chr_end=i'=>\$chr_end, 'image_type=s'=>\$image_type, 'force=s'=>\$force);

my $usage=<<INFO;

  Example usage: 
  
  CNView.pl  --reference_build=hg19  --cnv_file=/gscmnt/ams1184/info/model_data/2875816457/build111674790/variants/cnvs.hq  --working_dir=/gscmnt/sata132/techd/mgriffit/hg1/cnvs/  --sample_name=hg1
  
  Optional:
  --model_build               Instead of supplying a 'cnvs.hq' file, you can specify a somatic variation build ID and the file will be found automatically
  --gene_targets_file=PATH    List of gene symbols. These will be highlighted if they meet the minimum CNV cutoffs.  The cancer gene census will be used by default.
                              Example lists can be found here: /gscmnt/sata132/techd/mgriffit/reference_annotations/hg19/gene_symbol_lists/
  --name=NAME                 Name of gene symbol set to be highlighted ('Cancer Genes' will be used by default).
  --cnv_results_file=PATH     If you have already generated a CNV results file for a list of genes, you can supply it and skip to the graph generation step
  --ideogram_file=PATH        Chromosome ideogram coordinates for the reference build.  Will be automatically selected based on the reference build specified
                              An examples ideogram files can be found here: /gscmnt/sata132/techd/mgriffit/reference_annotations/hg19/ideogram/
  --chr                       If you wish to limit analysis to a single chromosome, use --chr (e.g. --chr=chr17). Otherwise all will be processed
  --chr_start, --chr_end      If you specify a single chromosome, you may also specify a range on that chromosome using --chr_start and --chr_end
  --image_type                Type of image files created. Default is 'jpeg'. Can also be 'png', 'bmp', 'tiff', 'pdf', or 'none' 
  --verbose=1                 More verbose output
  --force=1                   Force creation of subdirectories without prompting

INFO


#Set amplification / deletion cutoffs that will be used to created filtered files of results for convenience in downstream analyses
#Try two fold gain and 1/2 loss to start
my $amp_cutoff = 2;
my $del_cutoff = -0.5;


#Check/format all input parameters and options
#Based on the reference build specified, define paths to transcript/gene id mapping files
my ($gene_to_trans_file, $transcript_bed_file, $subdir, $outfile, $outfile_amp, $outfile_del, $outfile_ampdel, $name_f);
my $extra_flank = 100000;
&checkInput();

my $entrez_ensembl_data = &loadEntrezEnsemblData();

#If the user is supplying a pre-computed CNV file, the following steps will be skipped
my ($gt_map, $t_coords, $window_size, $cnvs, $targets);
if ($cnv_results_file){
  if($verbose){ print BLUE, "\n\nUsing user supplied CNV results file: $cnv_results_file", RESET; }
}else{
  if($verbose){ print BLUE, "\n\nGenerating CNV results file: $cnv_results_file", RESET; }

  #Load the gene to transcript name mappings.
  #Key on unique gene ID - reference to a hash containing all transcript IDs associated with that gene
  $gt_map = &loadGtMap();

  #Load the transcripts outer coords
  $t_coords = &loadTranscriptCoords();

  #Load the CNV window coordinates and copy number estimates
  $cnvs = &loadCnvData();

  #Load the gene targets of interest
  $targets = &loadTargetGenes();

  #For each gene target, get the corresponding transcripts
  #Merge these coordinates and get the grand outer coordinates for the gene of interest
  #Then get the mean CNV Diff for each gene
  &calculateMeanCnvDiff();

  #Print result to an output file
  &printCnvResultFile();
}

#Assume that the CNview.R script resides in the same dir as this script and obtain the path automatically

#Execute the R script that generates the CNV plots
my $rscript = "$script_dir/"."CNView.R";
my $r_cmd = "$rscript '$name' $cnv_file $cnv_results_file $ideogram_file $subdir $chr $chr_start $chr_end $image_type";
my $r_cmd_stdout = "$subdir"."CNView.R.stdout";
my $r_cmd_stderr = "$subdir"."CNView.R.stderr";
if ($verbose){
  print BLUE, "\n\nExecuting R code:\n$r_cmd\n", RESET;
}else{
  $r_cmd .= " 1>$r_cmd_stdout 2>$r_cmd_stderr";
}
Genome::Sys->shellcmd(cmd => $r_cmd);

#Annotate all gene entries in the output file with a cytoband label using the coordinates in the ideogram data file, and the coordinates of the gene

#Parse import the coordinates of the ideogram file using a subroutine
my $ideo_data = &importIdeogramData('-ideogram_file'=>$ideogram_file);

#Modify the output files ($cnv_results_file) above to annotate each gene with the corresponding cytoband(s) spanned by the gene
#Load file, get: (chr, start, end), determine cytobands that overlap, update output file, replace old file with new
#e.g. of cytoband: 11q13.2  OR  11q13.2 - 11q13.3
my @results_files = ($outfile, $outfile_amp, $outfile_del, $outfile_ampdel);
foreach my $results_file (@results_files){
  if (-e $results_file){
    #Parse input file and determine cytobands for each gene
    my %data;
    my $new_cnv_results_file = "$results_file".".tmp";
    my $header_line = '';
    my $header = 1;
    my $l = 0;
    open (CNV, "$results_file") || die "\n\nCould not open CNV results file: $results_file\n\n";
    my %columns;
    while(<CNV>){
      $l++;
      chomp($_);
      my @line = split("\t", $_);
      if ($header){
        $header = 0;
        $header_line = $_;
        my $p = 0;
        foreach my $col (@line){
          $columns{$col}{position} = $p;
          $p++;
        }
        next();
      }
      $data{$l}{line} = $_;
      my $chr = $line[$columns{'Chr'}{position}];
      my $chr_start = $line[$columns{'Start'}{position}];
      my $chr_end = $line[$columns{'End'}{position}];
      my $cytoband_string = &getCytoband('-ideo_data'=>$ideo_data, '-chr'=>$chr, '-chr_start'=>$chr_start, '-chr_end'=>$chr_end);
      $data{$l}{cytoband} = $cytoband_string;
    }
    close (CNV);

    #Print out new file containing the cytoband info
    open(CNV_NEW, ">$new_cnv_results_file") || die "\n\nCould not open new CNV results file: $new_cnv_results_file\n\n";
    print CNV_NEW "$header_line\tCytoband\n";
    foreach my $l (sort {$a <=> $b} keys %data){
      print CNV_NEW "$data{$l}{line}\t$data{$l}{cytoband}\n";
    }
    close(CNV_NEW);

    #Overwrite the old file with the new file
    my $mv_cmd = "mv $new_cnv_results_file $results_file";
    Genome::Sys->shellcmd(cmd => $mv_cmd);

  }else{
    print YELLOW, "\n\nCNV results file not found - not possible to annotate with Cytoband info\n\n", RESET;
  }
}

if ($verbose){ print BLUE, "\n\nResults written to:\n$subdir\n\n", RESET; }

exit();


#####################################################################################################################################################################
#Check/format all input parameters and options
#####################################################################################################################################################################
sub checkInput{
  
  #Check for required parameters
  unless ($reference_build && ($cnv_file || $model_build) && $working_dir && $sample_name){
    print GREEN, "$usage", RESET;
    print "\n\n";
    exit();
  }

  #Based on the reference build specified, define paths to transcript/gene id mapping files
  my $clinseq_annotations_dir = "/gscmnt/sata132/techd/mgriffit/reference_annotations/";
  if ($reference_build =~ /hg18/i){
    $gene_to_trans_file = $clinseq_annotations_dir . "hg18/transcript_to_gene/ALL.Genes.map";
    $transcript_bed_file = $clinseq_annotations_dir . "hg18/ALL.Genes.bed";
    $ideogram_file = $clinseq_annotations_dir . "hg18/ideogram/ChrBandIdeogram.tsv";
  }elsif($reference_build =~ /hg19/i){
    $gene_to_trans_file = $clinseq_annotations_dir . "hg19/transcript_to_gene/ALL.Genes.map";
    $transcript_bed_file = $clinseq_annotations_dir . "hg19/ALL.Genes.bed";
    $ideogram_file = $clinseq_annotations_dir . "hg19/ideogram/ChrBandIdeogram.tsv";
  }else{
    print RED, "\n\nSpecified reference build not understood.  Allowed values are: hg18, hg19\n\n", RESET;
    exit();
  }
  unless (-e $gene_to_trans_file && -e $transcript_bed_file && -e $ideogram_file){
    print RED, "\n\nOne or more of the following annotation files is missing:\n$gene_to_trans_file\n$transcript_bed_file\n$ideogram_file\n\n", RESET;
    exit();
  }
  #Set the default gene targets file if it wasnt specified by the user
  unless ($gene_targets_file){
    $gene_targets_file = $clinseq_annotations_dir . "GeneSymbolLists/CancerGeneCensusPlus_Sanger.txt"
  }
  unless (-e $gene_targets_file){
    print RED, "\n\nGene targets file not found:\n$gene_targets_file\n\n", RESET;
    exit();
  }
  #Check input CNV data
  if ($cnv_file){
    unless (-e $cnv_file){
      print RED, "\n\nCNV file missing:\n$cnv_file\n\n", RESET;
      exit();
    }
  }elsif ($model_build){
    #Use 'genome' code to find the directory for the somatic variation build using the model build ID
    #genome model build list --show data_directory 103590621
    my $genome_cmd = "genome model build list --show data_directory $model_build";
    my @result = `$genome_cmd`;
    chomp(@result);
    my $data_directory = '';
    foreach my $result (@result){
      if ($result =~ /^\//){
        $data_directory = $result;
      }
    }
    unless ($data_directory){
      print RED, "\n\nCould not find data_directory for build: $model_build using query\n$genome_cmd\n\n", RESET;
      exit();
    }
    #Possible names of the CNV data file.  copy_number_output.out  cno_copy_number.csv  cnvs.hq
    my $path1 = "$data_directory/"."copy_number_output.out";
    my $path2 = "$data_directory/"."cno_copy_number.csv";
    my $path3 = "$data_directory/"."cnvs.hq";
    my $path4 = "$data_directory/variants/"."copy_number_output.out";
    my $path5 = "$data_directory/variants/"."cno_copy_number.csv";
    my $path6 = "$data_directory/variants/"."cnvs.hq";

    my @path_list = ($path1, $path2, $path3, $path4, $path5, $path6);
    my $path_found = 0;
    foreach my $path (@path_list){
      if (-e $path){
        $cnv_file = $path;
        $path_found = 1;
      } 
    }
    unless($path_found){
      print RED, "\n\nCould not find a copy number results file in the directory: $data_directory\n\n", RESET;
      exit();
    }
  }

  unless (-e $working_dir && -d $working_dir){
    print YELLOW, "\n\nWorking dir not valid:\n$working_dir", RESET;

    if ($force){
      print YELLOW, "\n\nCreating dir: $working_dir\n\n", RESET;
      mkdir($working_dir);
    }else{
      #Ask the user if they want to create this dir
      print GREEN, "\n\nCreate this directory (y/n)? ", RESET;
      my $answer = '';
      $answer = <>;
      chomp($answer);
      if ($answer =~ /^y$|^yes$/i){
        print GREEN, "\n\nCreating dir: $working_dir\n\n", RESET;
        mkdir($working_dir);
      }else{
        print RED, "\n\nAborting then ...\n\n", RESET;
        exit();
      }
    }
  }
  unless ($working_dir =~ /\/$/){
    $working_dir .= "/";
  }
  #Get/Set the name of the comparison
  if ($name){
    $name_f = $name;
    $name_f =~ s/ //g;
  }else{
    $name = "Cancer Genes";
    $name_f = $name;
    $name_f =~ s/ //g;
  }
  #Create a subdirectory within the working directory
  $subdir = "$working_dir"."CNView_"."$name_f/";
  unless (-e $subdir){
    mkdir($subdir);
  }
 
  #Check the chromosome specified by the user for format (should be chrN).   If not specified, set $chr to 'ALL' and set $chr_start and $chr_end to 0
  chomp($chr);
  chomp($chr_start);
  chomp($chr_end);
  if ($chr_start || $chr_end){
    if (($chr_start =~ /^\d+$/) && ($chr_end =~ /^\d+$/)){
      unless ($chr_start <= $chr_end){
        print RED, "\n\n--chr_start ($chr_start) must be less than --chr_end ($chr_end)\n\n", RESET;
        exit();
      }
      $chr_start -= $extra_flank;
      $chr_end += $extra_flank;
      if ($chr_start < 1){
        $chr_start = 1;
      }
    }else{
      print RED, "\n\n--chr_start ($chr_start) and --chr_end ($chr_end) must both be integers", RESET;
      exit();
    }
    unless($chr){
      print RED, "\n\nIf you are going to define --chr_start and --chr_end, you must also define a chromosome with --chr\n\n", RESET;
      exit();
    }
  }else{
    $chr_start = 0;
    $chr_end = 0;
  }

  if ($chr){
    unless ($chr =~ /^chr/){
      $chr = "chr$chr";
    }
    if ($chr =~ /all/i){
      $chr = "ALL";
    }
  }else{
    $chr = "ALL";
    $chr_start = 0;
    $chr_end = 0;
  }

  #Name the output file
  if ($chr eq "ALL"){
    $outfile = "$subdir"."CNView_"."$name_f".".tsv";
    $outfile_amp = "$subdir"."CNView_"."$name_f".".amp.tsv";
    $outfile_del = "$subdir"."CNView_"."$name_f".".del.tsv";
    $outfile_ampdel = "$subdir"."CNView_"."$name_f".".ampdel.tsv";
  }elsif (($chr =~ /chr/i) && $chr_start && $chr_end){
    $outfile = "$subdir"."CNView_"."$name_f"."_"."$chr"."_"."$chr_start-$chr_end".".tsv";
    $outfile_amp = "$subdir"."CNView_"."$name_f"."_"."$chr"."_"."$chr_start-$chr_end".".amp.tsv";
    $outfile_del = "$subdir"."CNView_"."$name_f"."_"."$chr"."_"."$chr_start-$chr_end".".del.tsv";
    $outfile_ampdel = "$subdir"."CNView_"."$name_f"."_"."$chr"."_"."$chr_start-$chr_end".".ampdel.tsv";
  }else{
    $outfile = "$subdir"."CNView_"."$name_f"."_"."$chr".".tsv";
    $outfile_amp = "$subdir"."CNView_"."$name_f"."_"."$chr".".amp.tsv";
    $outfile_del = "$subdir"."CNView_"."$name_f"."_"."$chr".".del.tsv";
    $outfile_ampdel = "$subdir"."CNView_"."$name_f"."_"."$chr".".ampdel.tsv";
  }

  #If specified check the cnv results file
  if ($cnv_results_file){
    unless (-e $cnv_results_file){
      print RED, "\n\nCNV results file could not be found: $cnv_results_file\n\n", RESET;
      exit();
    }
  }
  #Check for valid image type.  If not specified, default to jpeg
  chomp($image_type);
  if ($image_type){
    unless ($image_type =~ /^jpeg$|^png$|^bmp$|^tiff$|^pdf$|^none$/){
      print RED, "\n\nimage_type: ($image_type) not recognized.  Must be one of: 'jpeg', 'png', 'bmp', 'tiff', 'pdf', 'none'\n\n", RESET;
      exit();
    }
  }else{
    $image_type = "jpeg";
  }

  return();
}

#####################################################################################################################################################################
#Load the gene to transcript name mappings.
#Key on unique gene ID - reference to a hash containing all transcript IDs associated with that gene
#####################################################################################################################################################################
sub loadGtMap{
  my %gt_map;
  open (GT, "$gene_to_trans_file") || die "\n\nCould not open input file: $gene_to_trans_file\n\n";
  while(<GT>){
    chomp($_);
    my @line = split("\t", $_);
    my $tid = $line[0];
    my $gstring = $line[1];
    unless ($gstring){
      next();
    }
    my @genes = split(",", $gstring);
    foreach my $gene (@genes){
      if ($gt_map{$gene}){ 
        my $trans_ref = $gt_map{$gene}{trans};
        $trans_ref->{$tid}=1;
      }else{
        my %trans;
        $trans{$tid}=1;
        $gt_map{$gene}{trans} = \%trans;
      }
    }
  }
  close(GT);
  return(\%gt_map);
}


#####################################################################################################################################################################
#Load the transcripts outer coords
#####################################################################################################################################################################
sub loadTranscriptCoords{
  my %t_coords;
  
  open (TC, "$transcript_bed_file") || die "\n\nCould not open input file: $gene_to_trans_file\n\n";
  while(<TC>){
    chomp($_);
    my @line = split("\t", $_);
    my $tid = $line[3];
    my $chr_input = $line[0];

    $t_coords{$tid}{chr} = $chr_input;
    $t_coords{$tid}{start} = $line[1];
    $t_coords{$tid}{end} = $line[2];
  }
  close(TC);

  return(\%t_coords);
}

#####################################################################################################################################################################
#Load the CNV window coordinates and copy number estimates
#####################################################################################################################################################################
sub loadCnvData{
  my %cnvs;
  open(CNV, "$cnv_file") || die "\n\nCould not open input file: $cnv_file\n\n";
  my $c = 0;
  my $p1 = 0;
  my $p2 = 0;
  while(<CNV>){
    chomp($_);
    if ($_ =~ /^\#|^CHR/){
      next();
    }
    $c++;
    my @line = split("\t", $_);
    my $chr_input = "chr"."$line[0]";
    my $start = $line[1];
    if ($c == 1){$p1 = $start;}
    if ($c == 2){$p2 = $start;}
    $cnvs{$chr_input}{$start}{tumor} = $line[2];
    $cnvs{$chr_input}{$start}{normal} = $line[3];
    $cnvs{$chr_input}{$start}{diff} = $line[4];
  }
  close(CNV);
  $window_size = $p2 - $p1;
  if ($verbose){
    print YELLOW, "\n\nDetected a CNV window size of $window_size bp.  Using this for overlap calculations", RESET;
  }
  return(\%cnvs);
}

#####################################################################################################################################################################
#Load the gene targets of interest
#####################################################################################################################################################################
sub loadTargetGenes{
  my %targets;
  open (TARGET ,"$gene_targets_file") || die "\n\nCould not open gene targets file\n\n";
  while(<TARGET>){
    chomp($_);
    if ($_ =~ /(\S+)/){
      #Skip "n/a"
      if ($1 eq "n/a"){
        next();
      }else{
        $targets{$1}{found} = 0;
      }
    }
  }
  close(TARGET);
  return(\%targets);
}


#####################################################################################################################################################################
#For each gene target, get the corresponding transcripts
#Merge these coordinates and get the grand outer coordinates for the gene of interest
#Then get the mean CNV Diff for each gene
#####################################################################################################################################################################
sub calculateMeanCnvDiff{

  foreach my $target (sort keys %{$targets}){
    if ($verbose){
      print BLUE, "\n\n$target", RESET;
    }
    $targets->{$target}->{mean_diff} = 0;
    if ($gt_map->{$target}){
      #NOTE: A target is only found if the NAME matches exactly.  This should be improved...
      #Get the transcript list of this gene
      my $trans_ref = $gt_map->{$target}->{trans};
      my $tcount = keys %{$trans_ref};
      if ($verbose){
        print BLUE, "\n\tFound $tcount transcripts", RESET;
      }

      #Merge the transcript coordinates.  Make sure all coordinates are from the same chromosome
      my %chr_list;
      my $grand_chr = '';
      my $grand_start = 1000000000000000000000000000000000000;
      my $grand_end = 0;

      #Check chromosomes of transcripts.  If there is more than one, chose the one with the most support
      foreach my $tid (keys %{$trans_ref}){
        my $chr = $t_coords->{$tid}->{chr};
        #Ignore transcripts from the haplotype chromosomes - Add to this list as they come up...
        if ($chr =~ /chr4_ctg9_hap1|chr6_apd_hap1|chr6_mann_hap4|chr6_qbl_hap6|chr6_ssto_hap7|chr6_mcf_hap5|chr6_cox_hap2|chr6_dbb_hap3|chr17_ctg5_hap1/){
          next();        
        }
        $chr_list{$chr}{count}++;
      }
      my $i = 0;
      foreach my $chr (sort {$chr_list{$b}{count} <=> $chr_list{$a}{count}} keys %chr_list){
        $i++;
        if ($i == 1){
          $grand_chr = $chr;
        }
      }

      #Now using the chr identified above, get the grand coords
      foreach my $tid (keys %{$trans_ref}){
        my $chr = $t_coords->{$tid}->{chr};
        unless ($chr eq $grand_chr){
          next();
        }
        my $start = $t_coords->{$tid}->{start};
        my $end = $t_coords->{$tid}->{end};
        if ($start < $grand_start){$grand_start = $start;}
        if ($end > $grand_end){$grand_end = $end;}
      }
      my $chr_count = keys %chr_list;
      my $grand_size = $grand_end - $grand_start;
      if ($chr_count == 1){
        if ($verbose){
          print BLUE, "\n\t$grand_chr:$grand_start-$grand_end ($grand_size)", RESET;
        }
      }else{
        my @chrs = keys %chr_list;
        if ($verbose){
          print YELLOW, "\n\tDid not find a single chromosome for this gene!  (@chrs) - chose the one with the most support (i.e. most transcripts)", RESET;
          print BLUE, "\n\t$grand_chr:$grand_start-$grand_end ($grand_size)", RESET;
          #print Dumper %chr_list;
        }
      }
      #Make sure CNV value exists for this chromosome
      unless ($cnvs->{$grand_chr}){
        if ($verbose){
          print YELLOW, "\n\tNo CNV data defined for this chromosome: $grand_chr", RESET;
        }
        next();
      }

      #Now identify the copy number windows that overlap this gene
      my %chr_cnvs = %{$cnvs->{$grand_chr}};
      my $window_count = keys %chr_cnvs;
      if ($verbose){
        print BLUE, "\n\tFound $window_count windows for this chromosome", RESET;
      }

      #Store values for later
      $targets->{$target}->{found} = 1;
      $targets->{$target}->{chr} = $grand_chr;
      $targets->{$target}->{start}= $grand_start;
      $targets->{$target}->{end} = $grand_end;
      $targets->{$target}->{size} = $grand_size;

      #Calculate the average copy number for the gene of interest
      my $overlaps = 0;
      my @diffs;
      foreach my $pos (sort {$a <=> $b} keys %chr_cnvs){
        my $cnv_start = $pos;
        my $cnv_end = $pos+$window_size;
        if (($cnv_start >= $grand_start && $cnv_start <= $grand_end) || ($cnv_end >= $grand_start && $cnv_end <= $grand_end) || ($cnv_start <= $grand_start && $cnv_end >= $grand_end)){
          $overlaps++;
          push (@diffs, $chr_cnvs{$pos}{diff});
        }
      }
      my $sum = 0;
      my $mean_diff = 0;
      foreach my $diff (@diffs){$sum+=$diff;}
      if ($overlaps > 0){
        $mean_diff = $sum/$overlaps;
      }
      if ($verbose){
        print BLUE, "\n\tFound $overlaps overlapping CNV windows with mean diff: $mean_diff", RESET;
      }
      $targets->{$target}->{mean_diff} = $mean_diff;
  
    }else{
      if ($verbose){
        print YELLOW, "\n\tCould not find any transcripts for this gene: $target", RESET;
      }
    }
  }
  return();
}


#####################################################################################################################################################################
#Print result to an output file
#####################################################################################################################################################################
sub printCnvResultFile{

  #Write four files, one with all results, one with amplifications passing a cutoff, one with deletions passing a cutoff, and one with both amplifications and deletions passing cutoffs
  open (OUT, ">$outfile") || die "\n\nCould not open outfile: $outfile for writting\n\n";
  open (OUT_AMP, ">$outfile_amp") || die "\n\nCould not open outfile: $outfile_amp for writting\n\n";
  open (OUT_DEL, ">$outfile_del") || die "\n\nCould not open outfile: $outfile_del for writting\n\n";
  open (OUT_AMPDEL, ">$outfile_ampdel") || die "\n\nCould not open outfile: $outfile_ampdel for writting\n\n";
  print OUT "Symbol\tmapped_gene_name\tSample\tChr\tStart\tEnd\tMean CNV Diff\n";
  print OUT_AMP "Symbol\tmapped_gene_name\tSample\tChr\tStart\tEnd\tMean CNV Diff\n";
  print OUT_DEL "Symbol\tmapped_gene_name\tSample\tChr\tStart\tEnd\tMean CNV Diff\n";
  print OUT_AMPDEL "Symbol\tmapped_gene_name\tSample\tChr\tStart\tEnd\tMean CNV Diff\n";

  foreach my $target (sort {abs($targets->{$b}->{mean_diff}) <=> abs($targets->{$a}->{mean_diff})} keys %{$targets}){
    if ($targets->{$target}->{found}){
      my $fixed_gene_name = &fixGeneName('-gene'=>$target, '-entrez_ensembl_data'=>$entrez_ensembl_data, '-verbose'=>$verbose);
      my $target_chr = $targets->{$target}->{chr};
      my $target_chr_start = $targets->{$target}->{start};
      my $target_chr_end = $targets->{$target}->{end};
      my $mean_diff = $targets->{$target}->{mean_diff};
      my $string = "$target\t$fixed_gene_name\t$sample_name\t$target_chr\t$target_chr_start\t$target_chr_end\t$mean_diff\n";
      #Name the output file
      if ($chr eq "ALL"){
        print OUT "$string";
        if ($mean_diff >= $amp_cutoff){
          print OUT_AMP "$string";
          print OUT_AMPDEL "$string";
        }
        if ($mean_diff <= $del_cutoff){
          print OUT_DEL "$string";
          print OUT_AMPDEL "$string";
        }
      }elsif (($chr =~ /chr/i) && $chr_start && $chr_end){
        unless (($target_chr eq $chr) && (($target_chr_start >= $chr_start && $target_chr_start <= $chr_end) || ($target_chr_end >= $chr_start && $target_chr_end <= $chr_end) || ($target_chr_start <= $chr_start && $target_chr_end >= $chr_end))){
          next();
        }
        print OUT "$string";
        if ($mean_diff >= $amp_cutoff){
          print OUT_AMP "$string";
          print OUT_AMPDEL "$string";
        }
        if ($mean_diff <= $del_cutoff){
          print OUT_DEL "$string";
          print OUT_AMPDEL "$string";
        }

      }else{
        unless ($target_chr eq $chr){
          next();
        }
        print OUT "$string";
        if ($mean_diff >= $amp_cutoff){
          print OUT_AMP "$string";
          print OUT_AMPDEL "$string";
        }
        if ($mean_diff <= $del_cutoff){
          print OUT_DEL "$string";
          print OUT_AMPDEL "$string";
        }
      }
    }
  }
  $cnv_results_file = $outfile;

  close(OUT);
  close(OUT_AMP);
  close(OUT_DEL);
  close(OUT_AMPDEL);

  return();
}




