#!/usr/bin/env genome-perl
#Written by Malachi Griffith
#For a group of ClinSeq models, get SVs and create a master table that merges all cases together
#Merge at the level of SV positions and then separately at the level of genes

use strict;
use warnings;
use Getopt::Long;
use Term::ANSIColor qw(:constants);
use Data::Dumper;
use above 'Genome';
use Genome::Model::ClinSeq::Util qw(:all);
use Genome::Model::ClinSeq::OriginalScripts::Converge qw(:all);

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
  
  convergeSvs.pl  --model_group_id='50714'  --outdir=/tmp/converge_svs/  --label='BRAF'  --verbose=1

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
my $positions_outfile = "$outdir"."$label"."_SVs_Merged_PositionLevel.tsv";
my $genes_outfile = "$outdir"."$label"."_SVs_Merged_GeneLevel.tsv";
my $positions_outfile_categorical = "$outdir"."$label"."_SVs_Merged_PositionLevel_Categorical.tsv";
my $genes_outfile_categorical = "$outdir"."$label"."_SVs_Merged_GeneLevel_Categorical.tsv";

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

my %svs;
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
  $model_list{$model_id}{final_name} = $final_name;

  #Find the appropriate SV file
  my $clinseq_sv_dir = $data_directory . "/" . $final_name . "/sv/";
  print BLUE, "\n$final_name\t$model_id\t$model_name\t$clinseq_sv_dir", RESET;

  my $sv_file = $clinseq_sv_dir . "CandidateSvCodingFusions.tsv";
  unless (-e $sv_file){
    print RED, "\n\nCould not find SV file: $sv_file\n\n", RESET;
    exit(1);
  }

  #Parse the SV file
  my $header = 1;
  my %columns;
  open (SV, "$sv_file") || die "\n\nCould not open SV file: $sv_file\n\n";
  while(<SV>){
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

    my $gene_pair	= $line[$columns{'gene_pair'}{pos}];
    my $gene1	= $line[$columns{'gene1'}{pos}];
    my $gene2	= $line[$columns{'gene2'}{pos}];
    my $coord1 = $line[$columns{'coord1'}{pos}];
    my $coord2 = $line[$columns{'coord2'}{pos}];
    my $coord = $coord1 . "_" . $coord2;
    my $mapped_gene_name1	= $line[$columns{'mapped_gene_name1'}{pos}];
    my $mapped_gene_name2	= $line[$columns{'mapped_gene_name2'}{pos}];
    my $pairoscope_tumor_reads = $line[$columns{'pairoscope_tumor_reads'}{pos}];
    my $pairoscope_normal_reads = $line[$columns{'pairoscope_normal_reads'}{pos}];

    #Merge to the level of distinct CTX SV positions...
    if ($svs{$coord}){
      $svs{$coord}{recurrence}++;
      my $cases_ref = $svs{$coord}{cases};
      $cases_ref->{$final_name}=1;
      $svs{$coord}{pairoscope_tumor_reads} += $pairoscope_tumor_reads;
      $svs{$coord}{pairoscope_normal_reads} += $pairoscope_normal_reads;
    }else{
      $svs{$coord}{recurrence} = 1;
      $svs{$coord}{gene_pair} = $gene_pair;
      $svs{$coord}{gene1} = $gene1;
      $svs{$coord}{gene2} = $gene2;
      $svs{$coord}{coord1} = $coord1;
      $svs{$coord}{coord2} = $coord2;
      $svs{$coord}{mapped_gene_name1} = $mapped_gene_name1;
      $svs{$coord}{mapped_gene_name2} = $mapped_gene_name2;
      $svs{$coord}{pairoscope_tumor_reads} = $pairoscope_tumor_reads;
      $svs{$coord}{pairoscope_normal_reads} = $pairoscope_normal_reads;
      $svs{$coord}{line} = $line;
      my %cases;
      $cases{$final_name}=1;
      $svs{$coord}{cases} = \%cases;
    }

    #Merge to the level of distinct gene name pairs
    if ($genes{$gene_pair}){
      $genes{$gene_pair}{total_mutation_count}++;
      my $positions_ref = $genes{$gene_pair}{positions};
      $positions_ref->{$coord}=1;
      my $cases_ref = $genes{$gene_pair}{cases};
      $cases_ref->{$final_name}=1;
    }else{
      $genes{$gene_pair}{mapped_gene_name1} = $mapped_gene_name1;
      $genes{$gene_pair}{mapped_gene_name2} = $mapped_gene_name2;
      $genes{$gene_pair}{total_mutation_count} = 1;
      my %positions;
      $positions{$coord} = 1;
      $genes{$gene_pair}{positions} = \%positions;
      my %cases;
      $cases{$final_name} = 1;
      $genes{$gene_pair}{cases} = \%cases;
    }
  }
  close(SV);
}
my $new_model_count = keys %model_list;
unless ($model_count == $new_model_count){
  print RED, "\n\nDiscrepancy between number of models in input model group and number of common/subject names found (redundant models?)\n\n", RESET;
  exit(1);
}


#Print the position level recurrence summary
open (OUT1, ">$positions_outfile") || die "\n\nCould not open output file for writing: $positions_outfile\n\n";
open (OUT2, ">$positions_outfile_categorical") || die "\n\nCould not open output file for writing: $positions_outfile_categorical\n\n";
print OUT1 "coord\trecurrence_count\tmutated_samples\tsum_pairoscope_tumor_reads\tsum_pairoscope_normal_reads\t$header_line\n";
print OUT2 "coord\tsample\n";
foreach my $coord (sort keys %svs){
  my $cases_ref = $svs{$coord}{cases};
  my @cases = keys %{$cases_ref};
  my @sort_cases = sort @cases;
  my $sort_cases_string = join(",", @sort_cases);
  print OUT1 "$coord\t$svs{$coord}{recurrence}\t$sort_cases_string\t$svs{$coord}{pairoscope_tumor_reads}\t$svs{$coord}{pairoscope_normal_reads}\t$svs{$coord}{line}\n";
  foreach my $case (@sort_cases){
    print OUT2 "$coord\t$case\n";
  }
}
close(OUT1);
close(OUT2);

#Print the gene level recurrence summary
open (OUT1, ">$genes_outfile") || die "\n\nCould not open output file for writing: $genes_outfile\n\n";
open (OUT2, ">$genes_outfile_categorical") || die "\n\nCould not open output file for writing: $genes_outfile_categorical\n\n";
print OUT1 "gene_pair\tmapped_gene_name1\tmapped_gene_name2\ttotal_mutation_count\tmutated_sample_count\tmutated_position_count\tmutated_samples\tmutated_positions\n";
print OUT2 "gene_pair\tsample\n";
foreach my $gene_pair (sort keys %genes){
  my $positions_ref = $genes{$gene_pair}{positions};
  my @positions = keys %{$positions_ref};
  my $positions_count = scalar(@positions);
  my @sort_positions = sort @positions;
  my $sort_positions_string = join (",", @sort_positions);
  my $cases_ref = $genes{$gene_pair}{cases};
  my @cases = keys %{$cases_ref};
  my @sort_cases = sort @cases;
  my $sort_cases_string = join(",", @sort_cases);
  my $cases_count = scalar(@cases);
  print OUT1 "$gene_pair\t$genes{$gene_pair}{mapped_gene_name1}\t$genes{$gene_pair}{mapped_gene_name2}\t$genes{$gene_pair}{total_mutation_count}\t$cases_count\t$positions_count\t$sort_cases_string\t$sort_positions_string\n";
  foreach my $case (@sort_cases){
    print OUT2 "$gene_pair\t$case\n";
  }
}
close(OUT1);
close(OUT2);

print BLUE, "\n\nWrote results to: $outdir\n\n", RESET;

exit();







