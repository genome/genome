#!/usr/bin/env genome-perl -w
#Written by Malachi Griffith

use DBI;
use strict;
use Data::Dumper;
use Getopt::Long;
use Benchmark;
use Term::ANSIColor qw(:constants);
use File::Basename;

#Required
my $model_groups = '';
my $model_ids = '';
my $results_dir = '';

GetOptions ('model_groups=s'=>\$model_groups, 'model_ids=s'=>\$model_ids, 'results_dir=s'=>\$results_dir);

my $usage=<<INFO;

  Example usage: 
  
  buildCufflinksFpkmMatrix.pl  --model_groups='11112,11113'  --results_dir=/gscmnt/sata132/techd/mgriffit/hg1/rna_seq/hg1_vs_brcs/

  buildCufflinksFpkmMatrix.pl  --model_ids='2877143142,2877143421,2877143448,2877143499,2877143501,2878871474,2878871444'  --results_dir=/gscmnt/sata132/techd/mgriffit/hg1/rna_seq/hg1_vs_brcs/
  
  Notes:
  For now this script assumes that Cufflinks was run with the '-G' option on the same genes GTF for each model specified...

  Details:
  --model_groups               Model group ID(s) for RNA-seq expression models (comma separate multiple IDs)
  --model_ids                  List of individual model IDs (comma separate multiple IDs)
  --results_dir                Where the resulting matrix files will be written

INFO

unless (($model_groups || $model_ids) && $results_dir){
  print GREEN, "$usage", RESET;
  exit();
}

if ($model_groups && $model_ids){
  print RED, "\n\nSpecify --model_groups OR --model_ids,  not both\n\n", RESET;
  exit();
}

#Get the model IDs associated with this model group
my %models;
if ($model_groups){
  my @model_groups = split(",", $model_groups);
  foreach my $model_group (@model_groups){
    my @result = `genome model-group member list --model-group-id $model_group --style text`;
    chomp(@result);
    foreach my $line (@result){
      if ($line =~ /^(\d+)/){
        $models{$1}{sample_name}="unknown";
      }
    }
  }
}else{
  my @model_ids = split(",", $model_ids);
  foreach my $model_id (@model_ids){
    $models{$model_id}{sample_name}="unknown";
  }
}

#For each model get the common name, sample name, last succeeded build dir, etc.
foreach my $model_id (sort {$a <=> $b} keys %models){
  print BLUE, "\nModel: $model_id", RESET;
  my $query = "genome model list --filter id=$model_id --show id,subject_name,individual_common_name,last_succeeded_build_directory  --noheaders  --style=csv";
  #print YELLOW, "\n\n$query", RESET;
  my $result = `$query`;
  chomp($result);
  my @line = split(",", $result);
  $models{$model_id}{sample_name} = $line[1];
  $models{$model_id}{individual_name} = $line[2];
  $models{$model_id}{build_dir} = $line[3];
  $models{$model_id}{process} = 0;

  #For each model, if a successful build dir is present, find the neccessary input files...
  if ($models{$model_id}{build_dir} =~ /NULL/i){
    print YELLOW, "\n\tBuild not complete for model: $model_id\t$models{$model_id}{individual_name}\t$models{$model_id}{sample_name}", RESET;
  }else{
    my $accepted_hits_bam = "$models{$model_id}{build_dir}/alignments/accepted_hits.bam";
    my $junctions_bed = "$models{$model_id}{build_dir}/alignments/junctions.bed";
    my $genes_fpkm = "$models{$model_id}{build_dir}/expression/genes.fpkm_tracking";
    my $transcripts_fpkm = "$models{$model_id}{build_dir}/expression/isoforms.fpkm_tracking";
    
    if (-e $accepted_hits_bam && -e $junctions_bed && -e $genes_fpkm && -e $transcripts_fpkm){
      print BLUE, "\n\tBuild looks good so far", RESET;
      $models{$model_id}{accepted_hits_bam} = $accepted_hits_bam;
      $models{$model_id}{junctions_bed} = $junctions_bed;
      $models{$model_id}{genes_fpkm} = $genes_fpkm;
      $models{$model_id}{transcripts_fpkm} = $transcripts_fpkm;
      $models{$model_id}{process} = 1;
    }else{
      $models{$model_id}{accepted_hits_bam} = '';
      $models{$model_id}{junctions_bed} = '';
      $models{$model_id}{genes_fpkm} = '';
      $models{$model_id}{transcripts_fpkm} = '';
      print YELLOW, "\n\tBuild is missing alignment or expression file:\n\t\t$accepted_hits_bam\n\t\t$junctions_bed\n\t\t$genes_fpkm\n\t\t$transcripts_fpkm", RESET;
    }
  }
}

#Build an FPKM matrix of the gene and transcript expression values
my $genes_fpkm_outfile = "genes_fpkm.tsv";
&buildMatrix('-file_key'=>"genes_fpkm", '-outfile'=>$genes_fpkm_outfile);

my $transcripts_fpkm_outfile = "transcripts_fpkm.tsv";
&buildMatrix('-file_key'=>"transcripts_fpkm", '-outfile'=>$transcripts_fpkm_outfile);

print "\n\n";

#Summarize the files and produce a sample cuffdiff command
my $bam_string = '';
my $individual_string = '';
foreach my $model_id (sort {$models{$a}->{individual_name} cmp $models{$b}->{individual_name}} keys %models){
  unless($models{$model_id}{process}){
    next();
  }
  print "\n$model_id\t$models{$model_id}{individual_name}\t$models{$model_id}{accepted_hits_bam}";
  $bam_string .= "$models{$model_id}{accepted_hits_bam} ";
  $individual_string .= "$models{$model_id}{individual_name},";
}

print "\n\nExample cuffdiff command:\n\ncuffdiff -o $results_dir/cuffdiff -M /gscmnt/ams1102/info/model_data/2772828715/build106409296/annotation_data/rRNA_MT_pseudogene.gtf -L $individual_string -p 4  /gscmnt/ams1102/info/model_data/2772828715/build106409296/annotation_data/all_sequences.gtf $bam_string\n\n";

#print Dumper %models;


exit();



######################################################################################################################
#Build an FPKM matrix of gene or transcript expression values                                                        #
######################################################################################################################
sub buildMatrix{
  my %args = @_;
  my $file_key = $args{'-file_key'};
  my $outfile = $args{'-outfile'};

  my %matrix;
  my %feature_anno;
  
  print BLUE, "\n\nBuilding matrix for: $file_key", RESET;
  my $data_col = "FPKM";
  my $gene_name_col = "gene_short_name";
  my $gene_id_col = "gene_id";
  my $locus_col = "locus";

  my $first_file = 1;
  my $ref_model_id;
  my $ref_model_count = 0;
  my $current_model_count = 0;
  foreach my $model_id (sort {$models{$a}->{individual_name} cmp $models{$b}->{individual_name}} keys %models){
    my $individual_name = $models{$model_id}{individual_name};
    unless ($models{$model_id}{process}){
      print YELLOW, "\n\tSkipping: $individual_name", RESET;
      next();
    }
    my $infile = $models{$model_id}{$file_key};
    print BLUE, "\n\tStoring: $individual_name ($infile)", RESET;

    open (EXP, $infile) || die "\n\nCould not open input file: $infile\n\n";
    my $header =1;
    my %columns;
    while (<EXP>){
      chomp($_);
      my @line = split("\t", $_);
      if ($header == 1){
        my $colpos = 0;
        foreach my $colname (@line){
          $columns{$colname}{colpos}=$colpos;
          $colpos++;
        }
        $header = 0;
        next();
      }
      my $feature_id = $line[0];

      unless ($columns{$data_col}){
        print RED, "\n\nData column: $data_col is not defined in the file: $infile\n\n", RESET;
        exit();
      }
      my $fpkm = $line[$columns{$data_col}{colpos}];

      #Use the first file as a reference.  Make sure all subsequent files contain the same number of ids and that the ids actually match
      if ($first_file){
        $ref_model_count++;
        $feature_anno{$feature_id}{gene_name} = $line[$columns{$gene_name_col}{colpos}];
        $feature_anno{$feature_id}{gene_id} = $line[$columns{$gene_id_col}{colpos}];
        $feature_anno{$feature_id}{locus} = $line[$columns{$locus_col}{colpos}];
      }else{
        $current_model_count++;
        unless (defined($matrix{$ref_model_id}{$feature_id})){
          print RED, "\n\nA feature was found in this file that was not in the reference file (first one considered)\n\n", RESET;
          exit();
        }
      }
      #Store the expression value in the matrix
      $matrix{$model_id}{$feature_id} = $fpkm;
    }
    close(EXP);
    unless ($first_file){
      unless ($ref_model_count == $current_model_count){
        print RED, "\n\nThe feature count for this file ($current_model_count) does not match that of the reference file ($ref_model_count)\n\n", RESET;
        exit();
      }
    }
    $first_file = 0;
    $ref_model_id = $model_id;
    $current_model_count = 0;

    my %current_hash = %{$matrix{$model_id}};
    my $features_stored = keys %current_hash;
    print BLUE, " - Stored $features_stored features...", RESET;
  }


  print BLUE, "\n\nWriting out matrix to $outfile", RESET;
  open (OUT, ">$outfile") || die "\n\nCould not open output file: $outfile\n\n";
  #Write out the header line
  my $header = "feature_id\t$gene_id_col\t$gene_name_col\t$locus_col";
  foreach my $model_id (sort {$models{$a}->{individual_name} cmp $models{$b}->{individual_name}} keys %models){
    unless ($models{$model_id}{process}){
      next();
    }
    my $individual_name = $models{$model_id}{individual_name};
    my $col_name = "$individual_name"."_"."$data_col";
    $header .= "\t$col_name";
  }
  print OUT "$header\n";

  #Go through each model and feature and write out the matrix
  my %feature_list = %{$matrix{$ref_model_id}};
  foreach my $feature_id (sort keys %feature_list){
    print OUT "$feature_id\t$feature_anno{$feature_id}{gene_id}\t$feature_anno{$feature_id}{gene_name}\t$feature_anno{$feature_id}{locus}";
    foreach my $model_id (sort {$models{$a}->{individual_name} cmp $models{$b}->{individual_name}} keys %models){
      unless ($models{$model_id}{process}){
        next();
      }
      print OUT "\t$matrix{$model_id}{$feature_id}";
    }
    print OUT "\n";
  }
  close(OUT);

  return();
}





