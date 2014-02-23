#!/usr/bin/perl
#Written by Malachi Griffith & Obi Griffith
#For a group of ClinSeq models, summarize the number of drug interactions found for the following event types:

#WGS SNVs, WGS Indels, WGS amplifications
#Exome SNVs, Exome Indels
#RNA-seq outlier genes

#Longer term, this concept should use the new DGIDB command line modules for getting interactions directly from the database

use strict;
use warnings;
use Getopt::Long;
use Term::ANSIColor qw(:constants);
use Data::Dumper;
use above 'Genome';
use Genome::Model::ClinSeq::OriginalScripts::Util qw(:all);
use Genome::Model::ClinSeq::Util qw(:all);
use Genome::Model::ClinSeq::OriginalScripts::Converge qw(:all);


my $build_ids = '';
my $model_ids = '';
my $model_group_id = '';
my $event_types_list = '';
my $dgidb_subdir_names = '';
my $filter_name = '';
my $outdir = '';
my $verbose = 0;

GetOptions ('build_ids=s'=>\$build_ids, 'model_ids=s'=>\$model_ids, 'model_group_id=s'=>\$model_group_id,
            'event_types_list=s'=>\$event_types_list, 'dgidb_subdir_names=s'=>\$dgidb_subdir_names, 'filter_name=s'=>\$filter_name, 
            'outdir=s'=>\$outdir, 'verbose=i'=>\$verbose);

my $usage=<<INFO;
  Example usage: 

  convergeDrugInteractionCounts.pl  --model_group_id='44083'  --event_types_list='all'  --dgidb_subdir_names='drugbank,santa_monica_lung'  --filter_name='default'  --outdir=/tmp  --verbose=1

  Specify *one* of the following as input (each model/build should be a ClinSeq model)
  --build_ids                  Comma separated list of specific build IDs
  --model_ids                  Comma separated list of specific model IDs
  --model_group_id             A single genome model group ID

  Additional parameters
  --dgidb_subdir_names         A comma-separated list of drug-gene database source subdir(s) (e.g. 'drugbank,santa_monica_lung')
  --filter_name                The name appended to each file indicating which filter was applied (e.g. 'default', 'antineo', etc.).
                               Note: If both drugbank and santa_monica_lung are specified a combination of antineo (for drugbank) and default (for santa_monica_lung) is used.

  --event_types_list           Specify a comma separated list of valid event types.  Use 'all' to include all types
                               Allowed event types: 'snv,indel,cnv_gain,rna_cufflinks_absolute,rna_tophat_absolute'

  --outdir                     Path of the a directory for output files
  --verbose                    More descriptive stdout messages

  Test Clinseq model groups:
  44083                        ClinSeq Model Group where all three data types are present (WGS, Exome, RNA-seq)

INFO

unless (($build_ids || $model_ids || $model_group_id) && $dgidb_subdir_names && $filter_name && $event_types_list && $outdir){
  print RED, "\n\nRequired parameter missing", RESET;
  print GREEN, "\n\n$usage", RESET;
  exit(1);
}

#Set output file names
$outdir = &checkDir('-dir'=>$outdir, '-clear'=>"no");
my $outfile = "$outdir/ClinSeq_WGS-EXOME-RNASEQ_drugabble_targets_summary.txt";
my $outfile2 = "$outdir/ClinSeq_WGS-EXOME-RNASEQ_drugabble_targets_summary_totals.txt";
my $outfile3 = "$outdir/ClinSeq_WGS-EXOME-RNASEQ_drugabble_targets_summary_genes.txt";
my $outfile4 = "$outdir/ClinSeq_WGS-EXOME-RNASEQ_drugabble_targets_summary_collapsed.txt";
open (SUMMARY, ">$outfile") or die "can't open $outfile for write\n";
open (SUMMARY2, ">$outfile4") or die "can't open $outfile4 for write\n";
open (TOTALS, ">$outfile2") or die "can't open $outfile2 for write\n";
open (GENES, ">$outfile3") or die "can't open $outfile3 for write\n";

#Check event types supplied by the user against the list supplied by the user ...
my @allowed_event_types =  qw (snv indel cnv_gain rna_cufflinks_absolute rna_tophat_absolute);
my @event_types;
if ($event_types_list =~ /all/i){
  @event_types = @allowed_event_types;
}else{
  my @supplied_event_types = split(",", $event_types_list);
  my %allowed_event_types;
  foreach my $et (@allowed_event_types){
    $et =~ s/^\s+|\s+$//g;
    $allowed_event_types{$et}=1;
  }
  foreach my $et (@supplied_event_types){
    $et =~ s/^\s+|\s+$//g;
    unless ($allowed_event_types{$et}){
      print RED, "\n\nEvent type ($et) is not in allowed event type list ( @allowed_event_types )\n\n", RESET;
      exit(1);
    }
    push(@event_types, $et);
  }
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

#Get files:
# - drug-gene interaction files for each event type
# - annotated files for each event type containing potentially druggable results
my $files = &getFiles('-models_builds'=>$models_builds, '-filter_name'=>$filter_name, '-dgidb_subdir_names'=>$dgidb_subdir_names);

#Known druggable genes
my $k_result = &parseKnownDruggableFiles('-files'=>$files, '-event_types'=>\@event_types);
my $k_g = $k_result->{'genes'}; #Known druggable genes
my $k_i = $k_result->{'interactions'}; #Known druggable gene interactions
my $data_type_sum = $k_result->{'data_type_sum'};
#print Dumper $data_type_sum;

#Print summary of druggable genes broken down by event and data type 
my %totals;
my %data_type_sum_collapsed;
print SUMMARY "patient\tevent_type\tdata_type\tgene\tdrug\tdgidb_source\n";
foreach my $patient (sort keys %{$data_type_sum}){
  foreach my $dgidb_source (sort keys %{$data_type_sum->{$patient}}){
    foreach my $event_type (sort keys %{$data_type_sum->{$patient}->{$dgidb_source}}){
      foreach my $data_type (sort keys %{$data_type_sum->{$patient}->{$dgidb_source}->{$event_type}}){
        foreach my $gene (sort keys %{$data_type_sum->{$patient}->{$dgidb_source}->{$event_type}->{$data_type}}){
          #Get drug list for gene
          my %drugs = %{$k_g->{$gene}->{drug_list}};
          foreach my $drug (sort keys %drugs){
            $totals{$patient}{$data_type}{$gene}++;
            print SUMMARY "$patient\t$event_type\t$data_type\t$gene\t$drug\t$dgidb_source\n";
            my $collapsed_string="$patient\t$event_type\t$data_type\t$gene\t$drug";
            $data_type_sum_collapsed{$collapsed_string}{$dgidb_source}=1;
          }
        }
      }
    }
  }
}

print SUMMARY2 "patient\tevent_type\tdata_type\tgene\tdrug\tdgidb_source\n";
foreach my $collapsed_string (sort keys %data_type_sum_collapsed){
  print SUMMARY2 "$collapsed_string\t",join(",",keys %{$data_type_sum_collapsed{$collapsed_string}}),"\n";
}

print TOTALS "patient\tdata_type\tdruggable_gene_count\n";
print GENES "patient\tdata_type\tdruggable_gene\n";
foreach my $patient (sort keys %totals){
  foreach my $data_type (sort keys %{$totals{$patient}}){
    my $gene_count = keys %{$totals{$patient}{$data_type}};
    print TOTALS "$patient\t$data_type\t$gene_count\n";
    foreach my $gene (sort keys %{$totals{$patient}{$data_type}}){
      print GENES "$patient\t$data_type\t$gene\n";
    }
  }
}

print "\n\n";
close SUMMARY;
close SUMMARY2;
close TOTALS;
close GENES;
exit();


############################################################################################################################
#Get input files to be parsed                                                                                              #
############################################################################################################################
sub getFiles{
  my %args = @_;
  my $models_builds = $args{'-models_builds'};
  my $filter_name = $args{'-filter_name'};
  #my $dgidb_subdir_name = $args{'-dgidb_subdir_name'};
  my @dgidb_subdir_names = split(",",$args{'-dgidb_subdir_names'});
  my %files;

  if ($verbose){print BLUE, "\n\nGet annotation files and drug-gene interaction files from these builds", RESET;}
  my %mb = %{$models_builds->{cases}};
  foreach my $dgidb_subdir_name (@dgidb_subdir_names){
    #Note: If both DrugBank and SantaMonica are specified as sources it perhaps does not make sense to use default for drug bank. 
    #SantaMonica is just cancer drugs so you would be merging apples and oranges.
    #Need to use the 'default' (i.e., unfiltered) setting for SantaMonica and the 'neoplastic' filter for drugbank.
    if ($dgidb_subdir_name eq 'drugbank'){$filter_name='antineo';}
    if ($dgidb_subdir_name eq 'santa_monica_lung'){$filter_name='default';}
    foreach my $c (keys %mb){
      my $b = $mb{$c}{build};
      my $m = $mb{$c}{model};
      my $build_directory = $b->data_directory;
      my $subject_name = $b->subject->name;
      my $subject_common_name = $b->subject->common_name;
      my $build_id = $b->id;

      #If the subject name is not defined, die
      unless ($subject_name){
        print RED, "\n\nCould not determine subject name for build: $build_id\n\n", RESET;
        exit(1);
      }

      my $final_name = "Unknown";
      if ($subject_name){$final_name = $subject_name;}
      if ($subject_common_name){$final_name = $subject_common_name;}
      if ($verbose){print BLUE, "\n\t$final_name\t$build_id\t$build_directory", RESET;}

      my $topdir = "$build_directory"."/$final_name/";

      #Some event types could have come from exome, wgs, or wgs_exome... depending on the event type allow these options and check in order

      #1.) Look for SNV files
      my @snv_subdir_options = qw (wgs exome);
      my $snv_drug_file_name = "snvs.hq.tier1.v1.annotated.compact."."$filter_name".".tsv";
      foreach my $dir_name (@snv_subdir_options){
        my $drug_file_path = $topdir . "snv/$dir_name/dgidb/$dgidb_subdir_name/$snv_drug_file_name";
        if (-e $drug_file_path){
          $files{$final_name}{$dgidb_subdir_name}{snv}{$dir_name}{drug_file_path} = $drug_file_path;
        }else{
          print RED, "\n\nCould not find SNV drug-gene file for $final_name ($subject_name - $subject_common_name) in:\n\t$build_directory\n\n", RESET;
          exit(1);
        }
      }

      #2.) Look for InDel files
      my @indel_subdir_options = qw (wgs exome);
      my $indel_drug_file_name = "indels.hq.tier1.v1.annotated.compact."."$filter_name".".tsv";
      foreach my $dir_name (@indel_subdir_options){
        my $drug_file_path = $topdir . "indel/$dir_name/dgidb/$dgidb_subdir_name/$indel_drug_file_name";
        #If both files are present, store for later
        if (-e $drug_file_path){
          $files{$final_name}{$dgidb_subdir_name}{indel}{$dir_name}{drug_file_path} = $drug_file_path;
        }else{
          print RED, "\n\nCould not find INDEL drug-gene file for $final_name ($subject_name - $subject_common_name) in:\n\t$build_directory\n\n", RESET;
          exit(1);
        }
      }

      #3.) Look for CNV gain files
      my $cnv_gain_drug_file_name = "cnv.AllGenes_Ensembl58.amp."."$filter_name".".tsv";
      my $drug_file_path1 = $topdir . "cnv/dgidb/$dgidb_subdir_name/$cnv_gain_drug_file_name";
      my $drug_file_path2 = $topdir . "cnv/cnview/dgidb/$dgidb_subdir_name/$cnv_gain_drug_file_name";

      if (-e $drug_file_path1){
        $files{$final_name}{$dgidb_subdir_name}{cnv_gain}{wgs}{drug_file_path} = $drug_file_path1;
      }elsif(-e $drug_file_path2){
        $files{$final_name}{$dgidb_subdir_name}{cnv_gain}{wgs}{drug_file_path} = $drug_file_path2;
      }else{
        print RED, "\n\nCould not find CNV drug-gene file for $final_name ($subject_name - $subject_common_name) in:\n\t$build_directory\n\n", RESET;
        exit(1);
      }

      #4.) Look for Cufflinks RNAseq outlier expression files 
      my $rna_cufflinks_drug_file_name = "isoforms.merged.fpkm.expsort.top1percent."."$filter_name".".tsv";
      my $drug_file_path = $topdir . "rnaseq/tumor/cufflinks_absolute/isoforms_merged/dgidb/$dgidb_subdir_name/$rna_cufflinks_drug_file_name";
      if (-e $drug_file_path){
        $files{$final_name}{$dgidb_subdir_name}{rna_cufflinks_absolute}{rnaseq}{drug_file_path} = $drug_file_path;
      }else{
        print RED, "\n\nCould not find Cufflinks Absolute drug-gene file for $final_name ($subject_name - $subject_common_name) in:\n\t$build_directory\n\t$drug_file_path\n\n", RESET;
        exit(1);
      }

      #5.) Look for Tophat junction RNAseq outlier expression files
      my $rna_tophat_drug_file_name = "Ensembl.Junction.GeneExpression.top1percent."."$filter_name".".tsv";
      $drug_file_path = $topdir . "rnaseq/tumor/tophat_junctions_absolute/dgidb/$dgidb_subdir_name/$rna_tophat_drug_file_name";
      if (-e $drug_file_path){
        $files{$final_name}{$dgidb_subdir_name}{rna_tophat_absolute}{rnaseq}{drug_file_path} = $drug_file_path;
      }else{
        print RED, "\n\nCould not find Tophat Absolute drug-gene file for $final_name ($subject_name - $subject_common_name) in:\n\t$build_directory\n\t$drug_file_path\n\n", RESET;
        exit(1);
      }
    }
  }
  return(\%files);
}

############################################################################################################################
#Parse known druggable files                                                                                               #
############################################################################################################################
sub parseKnownDruggableFiles{
  my %args = @_;
  my $files = $args{'-files'};
  my @event_types = @{$args{'-event_types'}};

  print BLUE, "\n\nParsing files containing known drug-gene interaction data", RESET;

	#Store all results organized by gene and separately by drug-gene interaction and separately by patient
  my %result;
  my %genes;
  my %interactions;
  my %patients;
  my %data_type_sum;

	#Note that the druggable files are already filtered down to only the variant affected genes with a drug interaction
  #To get a sense of the total number of events will have to wait until the annotation files are being proccessed 

  foreach my $patient (keys %{$files}){
    print BLUE, "\n\t$patient", RESET;
    foreach my $dgidb_source (keys %{$files->{$patient}}){
      foreach my $event_type (@event_types){
        foreach my $data_type (sort keys %{$files->{$patient}->{$dgidb_source}->{$event_type}}){
          my $drug_file_path = $files->{$patient}->{$dgidb_source}->{$event_type}->{$data_type}->{drug_file_path};
          print BLUE, "\n\t\t$drug_file_path", RESET;

          open (IN, "$drug_file_path") || die "\n\nCould not open gene-drug interaction file: $drug_file_path\n\n";
          my $header = 1;
          my %columns;
          while(<IN>){
            chomp($_);
            my @line = split("\t", $_);
            if ($header){
              my $p = 0;
              foreach my $column (@line){
                $columns{$column}{position} = $p;
                $p++;
              }
              $header = 0;
              next();
            }
            my $mapped_gene_name = $line[$columns{'mapped_gene_name'}{position}];
            my $drug_name = $line[$columns{'drug_name'}{position}];
            my $interaction = "$mapped_gene_name"."_"."$drug_name";
            my $drug_class = "unknown";
            if (defined($columns{'drug_class'})){
              $drug_class = $line[$columns{'drug_class'}{position}];
            }
            $interactions{$interaction}{mapped_gene_name} = $mapped_gene_name;
            $interactions{$interaction}{drug_name} = $drug_name;

            #Store drug class list from santa monica db
            if (defined($genes{$mapped_gene_name}{drug_class})){
              my $classes = $genes{$mapped_gene_name}{drug_class};
              $classes->{$drug_class} = 1;
            }else{
              my %classes;
              $classes{$drug_class} = 1;
              $genes{$mapped_gene_name}{drug_class} = \%classes;
            }

            #If the gene has any events it will be associated with all drugs that interact with that gene
            if (defined($genes{$mapped_gene_name}{drug_list})){
              my $drugs = $genes{$mapped_gene_name}{drug_list};
              $drugs->{$drug_name} = 1;
            }else{
              my %drugs;
              $drugs{$drug_name} = 1;
              $genes{$mapped_gene_name}{drug_list} = \%drugs;
            }

            #Store unique combinations of patient, event type, data type and gene (e.g., PNC4  snv  wgs  KRAS)
            $data_type_sum{$patient}{$dgidb_source}{$event_type}{$data_type}{$mapped_gene_name} = 1;

            #Add patient lists specific to this event type
            if (defined($genes{$mapped_gene_name}{$event_type})){
              my $patients = $genes{$mapped_gene_name}{$event_type}{patient_list};
              $patients->{$patient} = 1;
            }else{
              my %patients;
              $patients{$patient} = 1;
              $genes{$mapped_gene_name}{$event_type}{patient_list} = \%patients;
            }

            if (defined($interactions{$interaction}{$event_type})){
              my $patients = $interactions{$interaction}{$event_type}{patient_list};
              $patients->{$patient} = 1;
            }else{
              my %patients;
              $patients{$patient} = 1;
              $interactions{$interaction}{$event_type}{patient_list} = \%patients;
            }

            #Create or update the grand list of patients with ANY events hitting this gene
            if (defined($genes{$mapped_gene_name}{grand_list})){
              my $patients = $genes{$mapped_gene_name}{grand_list};
              $patients->{$patient} = 1;
            }else{
              my %patients;
              $patients{$patient} = 1;
              $genes{$mapped_gene_name}{grand_list} = \%patients;
            }

            if (defined($interactions{$interaction}{grand_list})){
              my $patients = $interactions{$interaction}{grand_list};
              $patients->{$patient} = 1;
            }else{
              my %patients;
              $patients{$patient} = 1;
              $interactions{$interaction}{grand_list} = \%patients;
            }
          }
          close(IN);
        }
      }
    }
  }
  $result{genes} = \%genes;
  $result{interactions} = \%interactions;
  $result{data_type_sum} = \%data_type_sum;
  return(\%result);
}


