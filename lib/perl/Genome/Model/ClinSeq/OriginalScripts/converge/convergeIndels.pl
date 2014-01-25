#!/usr/bin/env genome-perl
#Written by Malachi Griffith
#For a group of ClinSeq models, get Indels and create a master table that merges all cases together
#Merge at the level of INDEL positions and then separately at the level of genes

use strict;
use warnings;
use Getopt::Long;
use Term::ANSIColor qw(:constants);
use Data::Dumper;
use above 'Genome';
use Genome::Model::ClinSeq::Util qw(:all);
use Genome::Model::ClinSeq::OriginalScripts::Converge qw(:all);

#Required
my $build_ids = '';
my $model_ids = '';
my $model_group_id = '';
my $outdir = '';
my $label = '';
my $verbose = 0;
my $test = 0;

GetOptions ('build_ids=s'=>\$build_ids, 'model_ids=s'=>\$model_ids, 'model_group_id=s'=>\$model_group_id, 'outdir=s'=>\$outdir, 'label=s'=>\$label, 'verbose=i'=>\$verbose, 'test=i'=>\$test);

my $usage=<<INFO;
  Example usage: 
  
  convergeIndels.pl  --model_group_id='25307'  --outdir=/gscmnt/sata132/techd/mgriffit/braf_resistance/recurrence_indel_results/  --label='BRAF'  --verbose=1

  Specify *one* of the following as input (each model/build should be a ClinSeq model)
  --build_ids            Comma separated list of specific build IDs
  --model_ids            Comma separated list of specific model IDs
  --model_group_id       A single genome model group ID

  Combines INDEL results from a group of Clinseq models into a single report:
  --outdir               Path to directory for output files
  --label                Label to apply to output files as a prefix
  --verbose              More descriptive stdout messages
  --test                 Use --test=1 to limit BAM read counting to a small number of positions

INFO

unless (($build_ids || $model_ids || $model_group_id) && $outdir && $label){
  print RED, "\n\nRequired parameter missing", RESET;
  print GREEN, "\n\n$usage", RESET;
  exit(1);
}
unless ($outdir =~ /\/$/){
  $outdir .= "/";
}
unless (-e $outdir && -d $outdir){
  print RED, "\n\nOutput directory is not valid: $outdir\n\n", RESET;
  exit(1);
}

#Define paths/names of final output files
my $positions_outfile = "$outdir"."$label"."_INDELs_Merged_PositionLevel.tsv";
my $genes_outfile = "$outdir"."$label"."_INDELs_Merged_GeneLevel.tsv";
my $positions_outfile_categorical = "$outdir"."$label"."_INDELs_Merged_PositionLevel_Categorical.tsv";
my $genes_outfile_categorical = "$outdir"."$label"."_INDELs_Merged_GeneLevel_Categorical.tsv";

#Get the models/builds
my $models_builds;
if ($build_ids){
  my @build_ids = split(",", $build_ids);
  unless(scalar(@build_ids) > 0){
    print RED, "\n\nCould not parse build_ids list: $build_ids", RESET;
  }
  $models_builds = &getModelsBuilds('-builds'=>\@build_ids, '-verbose'=>$verbose);
}elsif($model_ids){
  my @model_ids = split(",", $model_ids);
  unless(scalar(@model_ids) > 0){
    print RED, "\n\nCould not parse model_ids list: $model_ids", RESET;
  }
  $models_builds = &getModelsBuilds('-models'=>\@model_ids, '-verbose'=>$verbose);
}elsif($model_group_id){
  $models_builds = &getModelsBuilds('-model_group_id'=>$model_group_id, '-verbose'=>$verbose);
}

my @models = @{$models_builds->{models}};
my @builds = @{$models_builds->{builds}};

my $model_count = scalar(@models);

my %indels;
my %genes;

#Cycle through the models and get their builds
my $header_line;
my %model_list;
foreach my $m (@models){
  my $model_name = $m->name;
  my $model_id = $m->genome_model_id;
  my $b = $m->last_succeeded_build || "NONE_FINISHED";
  unless ($b){
    print RED, "\n\nCould not find a succeeded build to use... for $model_name\n\n", RESET;
    next();
  }
  my $data_directory = $b->data_directory;
  my $patient = $m->subject;
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
  my $final_name;
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

  #Store model objects and values for later
  $model_list{$model_id}{model_name} = $model_name;
  $model_list{$model_id}{data_directory} = $data_directory;
  $model_list{$model_id}{patient} = $patient;
  $model_list{$model_id}{wgs_build} = $wgs_build;
  $model_list{$model_id}{exome_build} = $exome_build;
  $model_list{$model_id}{final_name} = $final_name;

  #Find the appropriate INDEL file
  my $clinseq_indel_dir = $data_directory . "/" . $final_name . "/indel/";
  if ($wgs_build && $exome_build){
    $clinseq_indel_dir .= "wgs_exome/";
  }elsif($wgs_build){
    $clinseq_indel_dir .= "wgs/";
  }elsif($exome_build){
    $clinseq_indel_dir .= "exome/";
  }
  print BLUE, "\n$final_name\t$model_id\t$model_name\t$clinseq_indel_dir", RESET;

  my $indel_file = $clinseq_indel_dir . "indels.hq.tier1.v1.annotated.compact.tsv";
  unless (-e $indel_file){
    print RED, "\n\nCould not find INDEL file: $indel_file\n\n", RESET;
    exit(1);
  }

  #Parse the INDEL file
  my $header = 1;
  my %columns;
  open (INDEL, "$indel_file") || die "\n\nCould not open INDEL file: $indel_file\n\n";
  while(<INDEL>){
    chomp($_);
    my $line = $_;
    my @line = split("\t", $line);
    if ($header == 1){
      $header_line = $line;
      $header = 0;
      my $p = 0;
      foreach my $head (@line){
        $columns{$head}{pos} = $p;
        $p++;
      }
      next();
    }

    my $coord = $line[$columns{'coord'}{pos}];
    my $gene_name = $line[$columns{'gene_name'}{pos}];
    my $mapped_gene_name = $line[$columns{'mapped_gene_name'}{pos}];
    my $ensembl_gene_id = $line[$columns{'ensembl_gene_id'}{pos}];
    my $aa_changes = $line[$columns{'aa_changes'}{pos}];
    my $ref_base = $line[$columns{'ref_base'}{pos}];
    my $var_base = $line[$columns{'var_base'}{pos}];
    my $gid = $ensembl_gene_id;

    #Merge to the level of distinct positions...
    if ($indels{$coord}){
      $indels{$coord}{recurrence}++;
      my $cases_ref = $indels{$coord}{cases};
      $cases_ref->{$final_name}=1;
    }else{
      $indels{$coord}{gene_name} = $gene_name;
      $indels{$coord}{mapped_gene_name} = $mapped_gene_name;
      $indels{$coord}{aa_changes} = $aa_changes;
      $indels{$coord}{ensembl_gene_id} = $ensembl_gene_id;
      $indels{$coord}{ref_base} = $ref_base;
      $indels{$coord}{var_base} = $var_base;
      $indels{$coord}{recurrence} = 1;
      $indels{$coord}{line} = $line;
      my %cases;
      $cases{$final_name}=1;
      $indels{$coord}{cases} = \%cases;
    }

    #Merge to the level of distinct gene names
    if ($genes{$gid}){
      $genes{$gid}{total_mutation_count}++;
      my $positions_ref = $genes{$gid}{positions};
      $positions_ref->{$coord}=1;
      my $cases_ref = $genes{$gid}{cases};
      $cases_ref->{$final_name}=1;
    }else{
      $genes{$gid}{gene_name} = $gene_name;
      $genes{$gid}{mapped_gene_name} = $mapped_gene_name;
      $genes{$gid}{ensembl_gene_id} = $ensembl_gene_id;
      $genes{$gid}{total_mutation_count} = 1;
      my %positions;
      $positions{$coord} = 1;
      $genes{$gid}{positions} = \%positions;
      my %cases;
      $cases{$final_name} = 1;
      $genes{$gid}{cases} = \%cases;
    }
  }
  close(INDEL);
}
my $new_model_count = keys %model_list;
unless ($model_count == $new_model_count){
  print RED, "\n\nDiscrepancy between number of models in input model group and number of common/subject names found (redundant models?)\n\n", RESET;
  exit();
}


#Print the position level recurrence summary
open (OUT1, ">$positions_outfile") || die "\n\nCould not open output file for writing: $positions_outfile\n\n";
open (OUT2, ">$positions_outfile_categorical") || die "\n\nCould not open output file for writing: $positions_outfile_categorical\n\n";
print OUT1 "coord\trecurrence_count\tmutated_samples\t$header_line\n";
print OUT2 "coord\tsample\n";
foreach my $coord (sort keys %indels){
  my $cases_ref = $indels{$coord}{cases};
  my @cases = keys %{$cases_ref};
  my @sort_cases = sort @cases;
  my $sort_cases_string = join(",", @sort_cases);
  print OUT1 "$coord\t$indels{$coord}{recurrence}\t$sort_cases_string\t$indels{$coord}{line}\n";
  foreach my $case (@sort_cases){
    print OUT2 "$coord\t$case\n";
  }
}
close(OUT1);
close(OUT2);

#Print the gene level recurrence summary
open (OUT1, ">$genes_outfile") || die "\n\nCould not open output file for writing: $genes_outfile\n\n";
open (OUT2, ">$genes_outfile_categorical") || die "\n\nCould not open output file for writing: $genes_outfile_categorical\n\n";
print OUT1 "ensembl_gene_id\tgene_name\tmapped_gene_name\ttotal_mutation_count\tmutated_sample_count\tmutated_position_count\tmutated_samples\tmutated_positions\n";
print OUT2 "ensembl_gene_id\tgene_name\tmapped_gene_name\tsample\n";
foreach my $gid (sort keys %genes){
  my $positions_ref = $genes{$gid}{positions};
  my @positions = keys %{$positions_ref};
  my $positions_count = scalar(@positions);
  my @sort_positions = sort @positions;
  my $sort_positions_string = join (",", @sort_positions);
  my $cases_ref = $genes{$gid}{cases};
  my @cases = keys %{$cases_ref};
  my @sort_cases = sort @cases;
  my $sort_cases_string = join(",", @sort_cases);
  my $cases_count = scalar(@cases);
  print OUT1 "$genes{$gid}{ensembl_gene_id}\t$genes{$gid}{gene_name}\t$genes{$gid}{mapped_gene_name}\t$genes{$gid}{total_mutation_count}\t$cases_count\t$positions_count\t$sort_cases_string\t$sort_positions_string\n";
  foreach my $case (@sort_cases){
    print OUT2 "$genes{$gid}{ensembl_gene_id}\t$genes{$gid}{gene_name}\t$genes{$gid}{mapped_gene_name}\t$case\n";
  }
}
close(OUT1);
close(OUT2);

#TODO:
#Calculate the WGS and/or Exome Variant Allele frequencies for all positions mutated in any sample, for all samples
#Use BAM read counts...


print BLUE, "\n\nWrote results to: $outdir\n\n", RESET;

exit();



