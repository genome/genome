#!/usr/bin/env genome-perl
#Written by Malachi Griffith
#For a group of ClinSeq models, CNV files.  Make it work for:
# - Somatic CNV differences summarized to the gene level)
# - Raw CNV values summarized at 10kb windows across the genome

#Summarize several aspects of the copy number results.  Produce a matrix for each of the following
#1.) Window level - simple matrix
#2.) Gene level - all genes
#3.) Gene level - amplified genes only (amplified in at least one build)
#4.) Gene level - deleted genes only

#Input:
#A list of Clinseq builds, models, or a Clinseq model-group

use strict;
use warnings;
use Getopt::Long;
use Term::ANSIColor qw(:constants);
use Data::Dumper;
use Statistics::Descriptive;
use above 'Genome';
use Genome::Model::ClinSeq::OriginalScripts::Util qw(:all);
use Genome::Model::ClinSeq::OriginalScripts::Converge qw(:all);

my $build_ids = '';
my $model_ids = '';
my $model_group_id = '';
my $amp_cutoff = '';
my $del_cutoff = '';
my $outdir = '';
my $verbose = 0;

GetOptions ('build_ids=s'=>\$build_ids, 'model_ids=s'=>\$model_ids, 'model_group_id=s'=>\$model_group_id, 
            'amp_cutoff=f'=>\$amp_cutoff, 'del_cutoff=f'=>\$del_cutoff,
            'outdir=s'=>\$outdir, 'verbose=i'=>\$verbose);


my $usage=<<INFO;
Example usage: 

  convergeStats.pl  --model_group_id='31779'  --outdir=/gscmnt/sata132/techd/mgriffit/luc/cnv/converge/  --verbose=1

  Specify *one* of the following as input (each model/build should be a ClinSeq model)
  --build_ids                Comma separated list of specific build IDs
  --model_ids                Comma separated list of specific model IDs
  --model_group_id           A singe genome model group ID

  --amp_cutoff               Min CNV Diff value that will be considered an amplification [default = 2]
  --del_cutoff               Min CNV Diff value that will be considered a deletion [default = -0.5]

  --outdir                   Path of the location to write output files to
  --verbose                  More descriptive stdout messages

  Test Clinseq model groups:
  31779                      ClinSeq - LUCs - v1

INFO

unless (($build_ids || $model_ids || $model_group_id) && $outdir){
  print RED, "\n\nRequired parameter missing", RESET;
  print GREEN, "\n\n$usage", RESET;
  exit(1);
}

#Get the models/builds
if ($verbose){print BLUE, "\n\nGet genome models/builds for supplied list", RESET;}
my $models_builds;
if ($build_ids){
    my @build_ids = split(",", $build_ids);
      $models_builds = &getModelsBuilds('-builds'=>\@build_ids, '-verbose'=>$verbose);
}elsif($model_ids){
    my @model_ids = split(",", $model_ids);
      $models_builds = &getModelsBuilds('-models'=>\@model_ids, '-verbose'=>$verbose);
}elsif($model_group_id){
    $models_builds = &getModelsBuilds('-model_group_id'=>$model_group_id, '-verbose'=>$verbose);
}else{
    print RED, "\n\nCould not obtains models/builds - check input to convergeCufflinksExpression.pl\n\n", RESET;
      exit();
}


#Check output dir
$outdir = &checkDir('-dir'=>$outdir, '-clear'=>"no");


#Set default gene level amplification/deletion cutoff if they were not specified by the user
unless ($amp_cutoff){
  $amp_cutoff = 2;
}
unless ($del_cutoff){
  $del_cutoff = -0.5;
}


my %files = %{&getCnvFiles('-models_builds'=>$models_builds, '-verbose'=>$verbose)};

#Get a list of column name that will be used
my %column_names;
foreach my $bid (sort {$files{$a}{column_name} cmp $files{$b}{column_name}} keys %files){
  my $column_name = $files{$bid}{column_name};
  $column_names{$column_name} = 1;
}


#Parse the raw CNV windows files and build a matrix of CNV value for all patients
my $cnvs;
$cnvs = &parseWindows('-files'=>\%files);

#Write new cnv files: cnvs.hq.diffs.matrix.tsv, cnvs.hq.median.tsv, cnvs.hq.mean.tsv
&writeCnvWindowOutputFiles('-cnvs'=>$cnvs, '-outdir'=>$outdir, '-column_names'=>\%column_names);

#Parse the gene level summarized values
my $genes = &parseGenes('-files'=>\%files);

&writeCnvGeneOutputFiles('-genes'=>$genes, '-outdir'=>$outdir, '-column_names'=>\%column_names, '-amp_cutoff'=>$amp_cutoff, '-del_cutoff'=>$del_cutoff);


if ($verbose){
  print "\n\nWrote output to: $outdir\n\n";
}

exit();



#######################################################################################################################
#Get the desired files from each model                                                                                #
#######################################################################################################################
sub getCnvFiles{
  my %args = @_;
  my $models_builds = $args{'-models_builds'};
  my $verbose = $args{'-verbose'};

  my %files;
  my $fc = 0;
  if ($verbose){print BLUE, "\n\nGet all Stats files within these builds that match 'cnvs.hq' AND 'cnv.AllGenes.tsv'", RESET;}
  my %mb = %{$models_builds->{cases}};
  foreach my $c (keys %mb){
    my $b = $mb{$c}{build};
    my $m = $mb{$c}{model};
    my $model_name = $m->name;
    my $model_id = $m->id;
    my $data_directory = $b->data_directory;
    my $subject_name = $b->subject->name;
    my $subject_common_name = $b->subject->common_name;
    my $build_id = $b->id;
    my $wgs_build = $b->wgs_build;
    my $exome_build = $b->exome_build;

    my ($wgs_common_name, $wgs_name, $exome_common_name, $exome_name);
    if ($wgs_build){
      $wgs_common_name = $wgs_build->subject->patient->common_name;
      $wgs_name = $wgs_build->subject->patient->name;
    }
    if ($exome_build){
      $exome_common_name = $exome_build->subject->patient->common_name;
      $exome_name = $exome_build->subject->patient->name;
    }
    #Get the patient common name from one of the builds, if none can be found, use the individual name instead, if that can't be found either set the name to 'UnknownName'
    my @names = ($wgs_common_name, $exome_common_name, $wgs_name, $exome_name);
    my $final_name = "Unknown";
    foreach my $name (@names){
      if ($name){
        $final_name = $name;
        last();
      }
    }
    
    unless ($final_name){
      print RED, "\n\nCould not determine a common or subject name from model: $model_name ($model_id)\n\n", RESET;
      exit();
    }

    if ($verbose){print BLUE, "\n\t$final_name\t$build_id\t$data_directory", RESET;}

    my @files;

    #Find the two input CNV data files in the build dir
    my $find_cmd1 = "find $data_directory/*/cnv/ -name cnvs.hq";
    if ($verbose){print YELLOW, "\n\t\t$find_cmd1", RESET;}
    my $cnvs_hq_file = `$find_cmd1`;
    chomp($cnvs_hq_file);

    my $cnvs_gene_file;

    my $find_cmd2 = "find $data_directory -name cnv.All_genes.tsv";
    if ($verbose){print YELLOW, "\n\t\t$find_cmd2", RESET;}
    my @tmp2 = `$find_cmd2`;
    chomp(@tmp2);
    foreach my $path (@tmp2){
      if ($path =~ /cnv\.All\_genes\.tsv/){
        $cnvs_gene_file = $path;
      }
    }

    if ($cnvs_hq_file && $cnvs_gene_file){
      $files{$build_id}{cnvs_hq_file} = $cnvs_hq_file;
      $files{$build_id}{cnvs_gene_file} = $cnvs_gene_file;
      $files{$build_id}{final_name} = $final_name;
      $files{$build_id}{subject_name} = $subject_name;
      $files{$build_id}{model_name} = $model_name;
    }else{
      if ($verbose){print RED, "\n\tCould not find cnvs.hq and  cnv.All_genes.tsv in build: $build_id ($data_directory)\n\n", RESET;}
      exit(1);
    }
  }

  #Determine a suitable name for data columns in the output file.  It must be distinct accross the builds.  Try these combinations in order:
  # final_name
  # final_name + subject_name
  # model_name
  # final_name + subject_name + build_id
  my $file_count =  keys %files;
  unless ($file_count > 1){
    print RED, "\n\nFound $file_count files to join (need at least 2 ...)\n\nLooked for files of name: Stats.tsv\n\n", RESET;
    exit(1);
  }
  my %finalnames;
  my %finalnames_subjectnames;
  my %modelnames;
  my %finalnames_subjectnames_builds;
  foreach my $bid (keys %files){
    my $final_name = $files{$bid}{final_name};
    my $subject_name = $files{$bid}{subject_name};
    my $model_name = $files{$bid}{model_name};
    my $id1 = $final_name;
    my $id2 = "$final_name"."_"."$subject_name";
    my $id3 = $model_name;
    my $id4 = "$final_name"."_"."$subject_name"."_"."$bid";
    
    $finalnames{$id1}=1;
    $finalnames_subjectnames{$id2}=1;
    $modelnames{$id3}=1;
    $finalnames_subjectnames_builds{$id4}=1;

    $files{$bid}{id1} = $id1;
    $files{$bid}{id2} = $id2;
    $files{$bid}{id3} = $id3;
    $files{$bid}{id4} = $id4;
  }

  my $id_choice;
  if ($file_count == keys %finalnames){
    $id_choice = "id1";
  }elsif($file_count == keys %finalnames_subjectnames){
    $id_choice = "id2";
  }elsif($file_count == keys %modelnames){
    $id_choice = "id3";
  }elsif($file_count == keys %finalnames_subjectnames_builds){
    $id_choice = "id4";
  }else{
    print RED, "\n\nCould not generate a suitable set of distinct labels for the data files to be joined...\n\n", RESET;
    exit();
  }

  foreach my $fc (keys %files){
    $files{$fc}{column_name} = $files{$fc}{"$id_choice"};
  }

  return(\%files);
}


######################################################################################################################
#Parse the raw CNV windows files and build a matrix of CNV value for all patients                                    #
######################################################################################################################
sub parseWindows{
  my %args = @_;
  my $files = $args{'-files'};

  if ($verbose){
    print BLUE, "\n\nImporting cnvs.hq data from each build", RESET;
  }

  my %cnvs;
  my $last_window_count = 0;
  foreach my $bid (sort {$files->{$a}->{column_name} cmp $files->{$b}->{column_name}} keys %{$files}){
    my $column_name = $files->{$bid}->{column_name};
    my $cnvs_hq_file = $files->{$bid}->{cnvs_hq_file};
    
    if ($verbose){
      print BLUE, "\n\tProcessing $column_name: $cnvs_hq_file", RESET;
    }
    
    my $window_count = 0;
    open (CNV, "$cnvs_hq_file") || die "\n\nCould not open input file: $cnvs_hq_file\n\n";
    my $header = 1;
    my %columns;
    my $o = 0;
    while(<CNV>){
      chomp($_);
      #Skip comment lines
      if ($_ =~ /^\#/){
        next();
      }
      my @line = split("\t", $_);
      if ($header){
        my $p = 0;
        foreach my $head (@line){
          $columns{$head}{p} = $p;
          $p++;
        }
        $header = 0;
        next();
      }
      $window_count++;
      $o++;
      my $chr = $line[$columns{'CHR'}{p}];
      my $pos = $line[$columns{'POS'}{p}];
      my $tumor = $line[$columns{'TUMOR'}{p}];
      my $normal = $line[$columns{'NORMAL'}{p}];
      my $diff = $line[$columns{'DIFF'}{p}];
      my $chr_pos = "$chr"."_"."$pos";
      $cnvs{$chr_pos}{order} = $o;
      $cnvs{$chr_pos}{chr} = $chr;
      $cnvs{$chr_pos}{pos} = $pos;
      $cnvs{$chr_pos}{$column_name}{diff} = $diff;
      $cnvs{$chr_pos}{$column_name}{tumor} = $tumor;
      $cnvs{$chr_pos}{$column_name}{normal} = $normal;

    }
    close(CNV);
    $files->{$bid}->{window_count} = $window_count;
    $last_window_count = $window_count;
  }

  #Make sure the window count was equal for all files
  my $stored_window_count = keys %cnvs;
  foreach my $bid (keys %{$files}){
    my $window_count = $files->{$bid}->{window_count};
    my $cnvs_hq_file = $files->{$bid}->{cnvs_hq_file};
    unless ($window_count == $stored_window_count){
      print RED, "\n\nNumber of windows stored ($stored_window_count) is not equal to window count for this file ($window_count): $cnvs_hq_file\n\n", RESET;
      exit(1);
    }
    unless ($window_count == $last_window_count){
      print RED, "\n\nNumber of windows is not equal across all cnvs.hq files in this group of builds!\n\n", RESET;
      exit(1);
    }
  }



  return(\%cnvs);
}


######################################################################################################################
#Parse the raw CNV windows files and build a matrix of CNV value for all patients                                    #
######################################################################################################################
sub parseGenes{
  my %args = @_;
  my $files = $args{'-files'};

  if ($verbose){
    print BLUE, "\n\nImporting cnvs gene data from each build", RESET;
  }

  my %genes;
  my $first_file = 1;
  foreach my $bid (sort {$files->{$a}->{column_name} cmp $files->{$b}->{column_name}} keys %{$files}){
    my $column_name = $files->{$bid}->{column_name};
    my $cnvs_gene_file = $files->{$bid}->{cnvs_gene_file};
    
    if ($verbose){
      print BLUE, "\n\tProcessing $column_name: $cnvs_gene_file", RESET;
    }

    my %columns;
    my $header = 1;
    my $o = 0;
    open (CNV, "$cnvs_gene_file") || die "\n\nCould not open CNV gene file: $cnvs_gene_file\n\n";
    while (<CNV>){
      chomp($_);
      my @line = split("\t", $_);
      if ($header){
        my $p = 0;
        foreach my $head_name (@line){
          $columns{$head_name}{pos} = $p;
          $p++;
        }
        $header = 0;
        next();
      }
      unless ($columns{'gene_id'} && $columns{'gene_name'} && $columns{'mapped_gene_name'} && $columns{'chr'} && $columns{'start'} && $columns{'end'} && $columns{'cytoband'} && $columns{'mean_cnv_diff'} && $columns{'cnvhmm_status'}){
        print RED, "\n\nCould not find a neccessary column in file: $cnvs_gene_file\n\n", RESET;
        exit(1);
      }

      #Grab data columns
      my $gid = $line[$columns{'gene_id'}{pos}];
      my $gene_name = $line[$columns{'gene_name'}{pos}];
      my $mapped_gene_name = $line[$columns{'mapped_gene_name'}{pos}];
      my $chr = $line[$columns{'chr'}{pos}];
      my $start = $line[$columns{'start'}{pos}];
      my $end = $line[$columns{'end'}{pos}];
      my $cytoband = $line[$columns{'cytoband'}{pos}];
      my $mean_cnv_diff = $line[$columns{'mean_cnv_diff'}{pos}];
      my $cnvhmm_status = $line[$columns{'cnvhmm_status'}{pos}];

      #If this is the first file processed, add this gene to the list.  If it is not, make sure the gene is already in the list and add the data only
      if ($first_file){
        $o++;
        $genes{$gid}{gene_name} = $gene_name;
        $genes{$gid}{mapped_gene_name} = $mapped_gene_name;
        $genes{$gid}{chr} = $chr;
        $genes{$gid}{start} = $start;
        $genes{$gid}{end} = $end;
        $genes{$gid}{cytoband} = $cytoband;
        $genes{$gid}{order} = $o;
        my %mean_cnv_diffs;
        $mean_cnv_diffs{$column_name}{mean_cnv_diff} = $mean_cnv_diff;
        $genes{$gid}{mean_cnv_diffs} = \%mean_cnv_diffs;

        my %cnvhmm_statuses;
        $cnvhmm_statuses{$column_name}{cnvhmm_status} = $cnvhmm_status;
        $genes{$gid}{cnvhmm_statuses} = \%cnvhmm_statuses;

      }else{
        unless (defined($genes{$gid})){
          print RED, "\n\nFound a gene gid ($gid) in this file that was not in the first file of this set\n\n", RESET;
          exit(1);
        }
        my $mean_cnv_diffs = $genes{$gid}{mean_cnv_diffs};
        $mean_cnv_diffs->{$column_name}->{mean_cnv_diff} = $mean_cnv_diff;
        my $cnvhmm_statuses = $genes{$gid}{cnvhmm_statuses};
        $cnvhmm_statuses->{$column_name}->{cnvhmm_status} = $cnvhmm_status;
      }
    }
    close (CNV);
    $first_file = 0;
  }
  return (\%genes);
}


######################################################################################################################
#Write new cnv files: cnvs.hq.diffs.matrix.tsv, cnvs.hq.median.tsv, cnvs.hq.mean.tsv                                 #
######################################################################################################################
sub writeCnvWindowOutputFiles{
  my %args = @_;
  my $cnvs = $args{'-cnvs'};
  my $outdir = $args{'-outdir'};
  my %column_names = %{$args{'-column_names'}};

  if ($verbose){
    print BLUE, "\n\nPrinting CNV window output files", RESET;
  }

  my $cnv_diff_matrix_file = $outdir . "cnvs.hq.diffs.matrix.tsv";
  my $cnv_hq_median_file = $outdir . "cnvs.hq.median.tsv";
  my $cnv_hq_mean_file = $outdir . "cnvs.hq.mean.tsv";

  my @columns = sort keys %column_names;
  my $columns_s = join("\t", @columns);
  my $header1 = "CHR\tPOS\t$columns_s";
  my $header2 = "CHR\tPOS\tTUMOR\tNORMAL\tDIFF";

  open (CNV1, ">$cnv_diff_matrix_file") || die "\n\nCould not open CNV diff matrix file for output: $cnv_diff_matrix_file\n\n";
  open (CNV2, ">$cnv_hq_median_file") || die "\n\nCould not open CNV diff medians file for output: $cnv_hq_median_file\n\n";
  open (CNV3, ">$cnv_hq_mean_file") || die "\n\nCould not open CNV diff means file for output: $cnv_hq_mean_file\n\n";
  print CNV1 "$header1\n";
  print CNV2 "$header2\n";
  print CNV3 "$header2\n";

  foreach my $chr_pos (sort {$cnvs->{$a}->{order} <=> $cnvs->{$b}->{order}} keys %{$cnvs}){
    my $chr = $cnvs->{$chr_pos}->{chr};
    my $pos = $cnvs->{$chr_pos}->{pos};
    my @tumor_readcounts;
    my @normal_readcounts;
    my @diffs;
    foreach my $column_name (sort keys %column_names){
      my $tumor = $cnvs->{$chr_pos}->{$column_name}->{tumor};
      my $normal = $cnvs->{$chr_pos}->{$column_name}->{normal};
      my $diff = $cnvs->{$chr_pos}->{$column_name}->{diff};
      push (@tumor_readcounts, $tumor);
      push (@normal_readcounts, $normal);
      push (@diffs, $diff);
    }
    my $diffs_s = join("\t", @diffs);
    print CNV1 "$chr\t$pos\t$diffs_s\n";

    #Calculate mean and median for: tumor read count, normal read count, and diff
    my $stat_tumor = Statistics::Descriptive::Full->new();
    $stat_tumor->add_data(@tumor_readcounts); 
    my $mean_tumor = $stat_tumor->mean();
    my $mean_tumor_p = sprintf("%.0f", $mean_tumor);
    my $median_tumor = $stat_tumor->median();
    my $median_tumor_p = sprintf("%.0f", $median_tumor);

    my $stat_normal = Statistics::Descriptive::Full->new();
    $stat_normal->add_data(@normal_readcounts); 
    my $mean_normal = $stat_normal->mean();
    my $mean_normal_p = sprintf("%.0f", $mean_normal);
    my $median_normal = $stat_normal->median();   
    my $median_normal_p = sprintf("%.0f", $median_normal);

    my $stat_diff = Statistics::Descriptive::Full->new();
    $stat_diff->add_data(@diffs); 
    my $mean_diff = $stat_diff->mean();
    my $mean_diff_p = sprintf("%.6f", $mean_diff);
    my $median_diff = $stat_diff->median();
    my $median_diff_p = sprintf("%.6f", $median_diff);

    print CNV2 "$chr\t$pos\t$median_tumor_p\t$median_normal_p\t$median_diff_p\n";
    print CNV3 "$chr\t$pos\t$mean_tumor_p\t$mean_normal_p\t$mean_diff_p\n";
  }

  close(CNV1);
  close(CNV2);
  close(CNV3);

  return();
}


######################################################################################################################
#Write new cnv files: cnv.genes.matrix.tsv, cnv.genes.amp.tsv, cnv.genes.del.tsv , cnv.genes.ampdel.tsv              #
######################################################################################################################
sub writeCnvGeneOutputFiles{
  my %args = @_;
  my $genes = $args{'-genes'};
  my $outdir = $args{'-outdir'};
  my %column_names = %{$args{'-column_names'}};
  my $amp_cutoff = $args{'-amp_cutoff'};
  my $del_cutoff = $args{'-del_cutoff'};

  if ($verbose){
    print BLUE, "\n\nPrint CNV gene output files", RESET;
  }

  my $cnv_genes_matrix_file = $outdir . "cnv.genes.matrix.tsv";
  my $cnv_genes_amp_file = $outdir . "cnv.genes.amp.tsv";
  my $cnv_genes_del_file = $outdir . "cnv.genes.del.tsv";
  my $cnv_genes_ampdel_file = $outdir . "cnv.genes.ampdel.tsv";

  my @columns = sort keys %column_names;
  my $columns_s = join("\t", @columns);

  foreach my $gid (sort {$genes->{$a}->{order} <=> $genes->{$b}->{order}} keys %{$genes}){

    my $mapped_gene_name = $genes->{$gid}->{mapped_gene_name};
    my $chr = $genes->{$gid}->{chr};
    my $start = $genes->{$gid}->{start};
    my $end = $genes->{$gid}->{end};
    my $cytoband = $genes->{$gid}->{cytoband};
    
    my @mean_cnv_diffs;
    my @amp_subjects;
    my @del_subjects;
    my @ampdel_subjects;

    foreach my $column_name (sort keys %column_names){
      my $mean_cnv_diffs = $genes->{$gid}->{mean_cnv_diffs};
      my $mean_cnv_diff = $mean_cnv_diffs->{$column_name}->{mean_cnv_diff};
      my $cnvhmm_statuses = $genes->{$gid}->{cnvhmm_statuses};
      my $cnvhmm_status = $cnvhmm_statuses->{$column_name}->{cnvhmm_status};

      push (@mean_cnv_diffs, $mean_cnv_diff);

      if ($mean_cnv_diff >= $amp_cutoff || $cnvhmm_status =~ /gain/i){
        push(@amp_subjects, $column_name);
        push(@ampdel_subjects, $column_name);
      }
      if ($mean_cnv_diff <= $del_cutoff || $cnvhmm_status =~ /loss/i){
        push(@del_subjects, $column_name);
        push(@ampdel_subjects, $column_name);
      }
    }
    my $mean_cnv_diffs_s = join("\t", @mean_cnv_diffs);
    $genes->{$gid}->{mean_cnv_diffs_s} = $mean_cnv_diffs_s;

    my $amp_subject_count = scalar(@amp_subjects);
    my $amp_subjects_s = "NA";
    if ($amp_subject_count > 0){
      $amp_subjects_s = join(",", @amp_subjects);
    }
    $genes->{$gid}->{amp_subject_count} = $amp_subject_count;
    $genes->{$gid}->{amp_subject_list} = $amp_subjects_s;

    my $del_subject_count = scalar(@del_subjects);
    my $del_subjects_s = "NA";
    if ($del_subject_count > 0){
      $del_subjects_s = join(",", @del_subjects);
    }
    $genes->{$gid}->{del_subject_count} = $del_subject_count;
    $genes->{$gid}->{del_subject_list} = $del_subjects_s;

    my $ampdel_subject_count = scalar(@ampdel_subjects);
    my $ampdel_subjects_s = "NA";
    if ($ampdel_subject_count > 0){
      $ampdel_subjects_s = join(",", @ampdel_subjects);
    }
    $genes->{$gid}->{ampdel_subject_count} = $ampdel_subject_count;
    $genes->{$gid}->{ampdel_subject_list} = $ampdel_subjects_s;
  }

  open (CNV1, ">$cnv_genes_matrix_file") || die "\n\nCould not open CNV gene diff file for output: $cnv_genes_matrix_file\n\n";
  open (CNV2, ">$cnv_genes_amp_file")    || die "\n\nCould not open CNV gene diff file for output: $cnv_genes_amp_file\n\n";
  open (CNV3, ">$cnv_genes_del_file")    || die "\n\nCould not open CNV gene diff file for output: $cnv_genes_del_file\n\n";
  open (CNV4, ">$cnv_genes_ampdel_file") || die "\n\nCould not open CNV gene diff file for output: $cnv_genes_ampdel_file\n\n";
  
  my $header = "ensembl_gene_id\tgene_name\tmapped_gene_name\tchr\tstart\tend\tcytoband";
  print CNV1 "$header\t$columns_s\n";
  print CNV2 "$header\tAmplificationSubjectCount\tAmplificationSubjectList\t$columns_s\n";
  print CNV3 "$header\tDeletionSubjectCount\tDeletionSubjectList\t$columns_s\n";
  print CNV4 "$header\tAmpDelSubjectCount\tAmpDelSubjectList\t$columns_s\n";

  foreach my $gid (sort {$genes->{$a}->{order} <=> $genes->{$b}->{order}} keys %{$genes}){
    print CNV1 "$gid\t$genes->{$gid}->{gene_name}\t$genes->{$gid}->{mapped_gene_name}\t$genes->{$gid}->{chr}\t$genes->{$gid}->{start}\t$genes->{$gid}->{end}\t$genes->{$gid}->{cytoband}\t$genes->{$gid}->{mean_cnv_diffs_s}\n";
  }
  foreach my $gid (sort {$genes->{$b}->{amp_subject_count} <=> $genes->{$a}->{amp_subject_count}} keys %{$genes}){
    if ($genes->{$gid}->{amp_subject_count} > 0){
    print CNV2 "$gid\t$genes->{$gid}->{gene_name}\t$genes->{$gid}->{mapped_gene_name}\t$genes->{$gid}->{chr}\t$genes->{$gid}->{start}\t$genes->{$gid}->{end}\t$genes->{$gid}->{cytoband}\t$genes->{$gid}->{amp_subject_count}\t$genes->{$gid}->{amp_subject_list}\t$genes->{$gid}->{mean_cnv_diffs_s}\n";
    }
  }
  foreach my $gid (sort {$genes->{$b}->{del_subject_count} <=> $genes->{$a}->{del_subject_count}} keys %{$genes}){
    if ($genes->{$gid}->{del_subject_count} > 0){
    print CNV3 "$gid\t$genes->{$gid}->{gene_name}\t$genes->{$gid}->{mapped_gene_name}\t$genes->{$gid}->{chr}\t$genes->{$gid}->{start}\t$genes->{$gid}->{end}\t$genes->{$gid}->{cytoband}\t$genes->{$gid}->{del_subject_count}\t$genes->{$gid}->{del_subject_list}\t$genes->{$gid}->{mean_cnv_diffs_s}\n";
    }
  }
  foreach my $gid (sort {$genes->{$b}->{ampdel_subject_count} <=> $genes->{$a}->{ampdel_subject_count}} keys %{$genes}){
    if ($genes->{$gid}->{ampdel_subject_count} > 0){
    print CNV4 "$gid\t$genes->{$gid}->{gene_name}\t$genes->{$gid}->{mapped_gene_name}\t$genes->{$gid}->{chr}\t$genes->{$gid}->{start}\t$genes->{$gid}->{end}\t$genes->{$gid}->{cytoband}\t$genes->{$gid}->{ampdel_subject_count}\t$genes->{$gid}->{ampdel_subject_list}\t$genes->{$gid}->{mean_cnv_diffs_s}\n";
    }
  }

  close(CNV1);
  close(CNV2);
  close(CNV3);
  close(CNV4);

  return();
}

