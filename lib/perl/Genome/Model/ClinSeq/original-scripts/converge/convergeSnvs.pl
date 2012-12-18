#!/usr/bin/env genome-perl
#Written by Malachi Griffith
#For a group of ClinSeq models, get SNVs and create a master table that merges all cases together
#Merge at the level of SNV positions and then separately at the level of genes

use strict;
use warnings;
use Getopt::Long;
use Term::ANSIColor qw(:constants);
use Data::Dumper;
use above 'Genome';
use Genome::Model::ClinSeq::Util qw(:all);
use Genome::Model::ClinSeq::Converge qw(:all);

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
  
  convergeSnvs.pl  --model_group_id='25307'  --outdir=/gscmnt/sata132/techd/mgriffit/braf_resistance/recurrence_snv_results/  --label='BRAF'  --verbose=1

  Specify *one* of the following as input (each model/build should be a ClinSeq model)
  --build_ids            Comma separated list of specific build IDs
  --model_ids            Comma separated list of specific model IDs
  --model_group_id       A single genome model group ID

  Combines SNV results from a group of Clinseq models into a single report:
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

#Hardcode Ensembl version parameter (should not matter since we are not using RNA-seq data here) - but must be defined for bam read counts code I'm using
my $ensembl_version = 58;

#Create a sub-directory for BAM read count results
my $bam_rc_subdir = $outdir . "bam_rc/";
mkdir($bam_rc_subdir);

#Define paths/names of final output files
my $positions_outfile = "$outdir"."$label"."_SNVs_Merged_PositionLevel.tsv";
my $genes_outfile = "$outdir"."$label"."_SNVs_Merged_GeneLevel.tsv";
my $positions_outfile_categorical = "$outdir"."$label"."_SNVs_Merged_PositionLevel_Categorical.tsv";
my $genes_outfile_categorical = "$outdir"."$label"."_SNVs_Merged_GeneLevel_Categorical.tsv";
my $positions_list = "$outdir"."$label"."_Master_SNV_List.tsv";
my $vaf_file_matrix_wgs = "$outdir"."$label"."_WGS_SNV_VAFs_Matrix.tsv";
my $mutation_status_file_matrix_wgs = "$outdir"."$label"."_WGS_SNV_MutationStatus_Matrix.tsv";
my $vaf_file_matrix_exome = "$outdir"."$label"."_Exome_SNV_VAFs_Matrix.tsv";
my $mutation_status_file_matrix_exome = "$outdir"."$label"."_Exome_SNV_MutationStatus_Matrix.tsv";

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

my %snvs;
my %genes;

#Cycle through the models and get their builds
my $header_line;
my %model_list;
my %wgs_sample_list;
my %exome_sample_list;
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
  if ($wgs_build){
    $wgs_sample_list{$final_name} = 1;
  }
  if ($exome_build){
    $exome_sample_list{$final_name} = 1;
  }

  #Store model objects and values for later
  $model_list{$model_id}{model_name} = $model_name;
  $model_list{$model_id}{data_directory} = $data_directory;
  $model_list{$model_id}{patient} = $patient;
  $model_list{$model_id}{wgs_build} = $wgs_build;
  $model_list{$model_id}{exome_build} = $exome_build;
  $model_list{$model_id}{final_name} = $final_name;

  #Find the appropriate SNV file
  my $clinseq_snv_dir = $data_directory . "/" . $final_name . "/snv/";
  if ($wgs_build && $exome_build){
    $clinseq_snv_dir .= "wgs_exome/";
  }elsif($wgs_build){
    $clinseq_snv_dir .= "wgs/";
  }elsif($exome_build){
    $clinseq_snv_dir .= "exome/";
  }
  print BLUE, "\n$final_name\t$model_id\t$model_name\t$clinseq_snv_dir", RESET;

  my $snv_file = $clinseq_snv_dir . "snvs.hq.tier1.v1.annotated.compact.tsv";
  unless (-e $snv_file){
    print RED, "\n\nCould not find SNV file: $snv_file\n\n", RESET;
    exit(1);
  }

  #Parse the SNV file
  my $header = 1;
  my %columns;
  open (SNV, "$snv_file") || die "\n\nCould not open SNV file: $snv_file\n\n";
  while(<SNV>){
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

    #Merge to the level of distinct positions...
    if ($snvs{$coord}){
      $snvs{$coord}{recurrence}++;
      my $cases_ref = $snvs{$coord}{cases};
      $cases_ref->{$final_name}=1;
    }else{
      $snvs{$coord}{gene_name} = $gene_name;
      $snvs{$coord}{mapped_gene_name} = $mapped_gene_name;
      $snvs{$coord}{ensembl_gene_id} = $ensembl_gene_id;
      $snvs{$coord}{aa_changes} = $aa_changes;
      $snvs{$coord}{ref_base} = $ref_base;
      $snvs{$coord}{var_base} = $var_base;
      $snvs{$coord}{recurrence} = 1;
      $snvs{$coord}{line} = $line;
      my %cases;
      $cases{$final_name}=1;
      $snvs{$coord}{cases} = \%cases;
    }

    #Merge to the level of distinct gene names
    if ($genes{$gene_name}){
      $genes{$gene_name}{total_mutation_count}++;
      my $positions_ref = $genes{$gene_name}{positions};
      $positions_ref->{$coord}=1;
      my $cases_ref = $genes{$gene_name}{cases};
      $cases_ref->{$final_name}=1;
    }else{
      $genes{$gene_name}{mapped_gene_name} = $mapped_gene_name;
      $genes{$gene_name}{total_mutation_count} = 1;
      my %positions;
      $positions{$coord} = 1;
      $genes{$gene_name}{positions} = \%positions;
      my %cases;
      $cases{$final_name} = 1;
      $genes{$gene_name}{cases} = \%cases;
    }
  }
  close(SNV);
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
foreach my $coord (sort keys %snvs){
  my $cases_ref = $snvs{$coord}{cases};
  my @cases = keys %{$cases_ref};
  my @sort_cases = sort @cases;
  my $sort_cases_string = join(",", @sort_cases);
  print OUT1 "$coord\t$snvs{$coord}{recurrence}\t$sort_cases_string\t$snvs{$coord}{line}\n";
  foreach my $case (@sort_cases){
    print OUT2 "$coord\t$case\n";
  }
}
close(OUT1);
close(OUT2);

#Print the gene level recurrence summary
open (OUT1, ">$genes_outfile") || die "\n\nCould not open output file for writing: $genes_outfile\n\n";
open (OUT2, ">$genes_outfile_categorical") || die "\n\nCould not open output file for writing: $genes_outfile_categorical\n\n";
print OUT1 "gene_name\tmapped_gene_name\ttotal_mutation_count\tmutated_sample_count\tmutated_position_count\tmutated_samples\tmutated_positions\n";
print OUT2 "gene_name\tsample\n";
foreach my $gene_name (sort keys %genes){
  my $positions_ref = $genes{$gene_name}{positions};
  my @positions = keys %{$positions_ref};
  my $positions_count = scalar(@positions);
  my @sort_positions = sort @positions;
  my $sort_positions_string = join (",", @sort_positions);
  my $cases_ref = $genes{$gene_name}{cases};
  my @cases = keys %{$cases_ref};
  my @sort_cases = sort @cases;
  my $sort_cases_string = join(",", @sort_cases);
  my $cases_count = scalar(@cases);
  print OUT1 "$gene_name\t$genes{$gene_name}{mapped_gene_name}\t$genes{$gene_name}{total_mutation_count}\t$cases_count\t$positions_count\t$sort_cases_string\t$sort_positions_string\n";
  foreach my $case (@sort_cases){
    print OUT2 "$gene_name\t$case\n";
  }
}
close(OUT1);
close(OUT2);

#Calculate the WGS and/or Exome Variant Allele frequencies for all positions mutated in any sample, for all samples


#First build a consolidated SNV positions file of the format:
#5:112176318-112176318   APC     APC     p.R1676T	    G	    C
open (POS, ">$positions_list") || die "\n\nCould not open master positions list for output: $positions_list\n\n";
print POS "coord\tgene_name\tmapped_gene_name\tensembl_gene_id\taa_changes\tref_base\tvar_base\n";
my $c = 0;
foreach my $coord (sort keys %snvs){
  $c++;
  print POS "$coord\t$snvs{$coord}{gene_name}\t$snvs{$coord}{mapped_gene_name}\t$snvs{$coord}{ensembl_gene_id}\t$snvs{$coord}{aa_changes}\t$snvs{$coord}{ref_base}\t$snvs{$coord}{var_base}\n";
  if ($test){
    if ($c >= 10){
      last();
    }
  }
}
close(POS);

#Now cycle through the build in the model group and run BAM read count jobs
#Cycle through the models and get their builds
print BLUE, "\n\nGetting BAM read counts for the master list of SNV positions for both WGS and Exome from all models in the model-group", RESET;
foreach my $model_id (sort keys %model_list){
  my $model_name = $model_list{$model_id}{model_name};
  my $data_directory = $model_list{$model_id}{data_directory};
  my $patient = $model_list{$model_id}{patient};
  my $wgs_build = $model_list{$model_id}{wgs_build};
  my $exome_build = $model_list{$model_id}{exome_build};
  my $final_name = $model_list{$model_id}{final_name};
  my $output_file = $bam_rc_subdir . "$final_name"."_BamReadCounts.tsv";
  $model_list{$model_id}{read_counts_file} = $output_file; 

  my $read_counts_script = "genome-perl `which genome` model clin-seq get-bam-read-counts";
  my $bam_rc_cmd = "$read_counts_script  --positions-file=$positions_list  --ensembl-version=$ensembl_version  --output-file=$output_file  --verbose=$verbose";
  $bam_rc_cmd .= "  --wgs-som-var-build=" . $wgs_build->id if $wgs_build;
  $bam_rc_cmd .= "  --exome-som-var-build=" . $exome_build->id if $exome_build;
  unless ($verbose){
    $bam_rc_cmd .= "  1>/dev/null 2>/dev/null";
  }
  print BLUE, "\n\t$bam_rc_cmd", RESET;
  Genome::Sys->shellcmd(cmd => $bam_rc_cmd);
}

#Now parse the read counts files and build a hash of SNVs and their variant allele frequencies (tumor and exome) for each subject
my %rc;
foreach my $model_id (sort keys %model_list){
  my $rc_file = $model_list{$model_id}{read_counts_file}; 
  my $wgs_build = $model_list{$model_id}{wgs_build};
  my $exome_build = $model_list{$model_id}{exome_build};
  my $final_name = $model_list{$model_id}{final_name};
  my $header = 1;
  my %columns;
  open (RC, "$rc_file") || die "\n\nCould not open read counts file: $rc_file\n\n";
  while(<RC>){
    chomp($_);
    my @line = split("\t", $_);
    if ($header){
      my $p = 0;
      foreach my $col (@line){
        $columns{$col}{position} = $p;
        $p++;
      }
      $header = 0;
      next();
    }
    my $coord = $line[$columns{'coord'}{position}];
    if ($wgs_build){
      $rc{$coord}{$final_name}{wgs_tumor_vaf} = $line[$columns{'WGS_Tumor_VAF'}{position}];
    }
    if ($exome_build){
      $rc{$coord}{$final_name}{exome_tumor_vaf} = $line[$columns{'Exome_Tumor_VAF'}{position}];
    }
  }
  close(RC);
}

#Write out the matrix files for both VAFs and mutation status (VAF > 0) in each subject at each position
#WGS VAFs
my @wgs_sample_list = sort keys %wgs_sample_list;
my $wgs_sample_count = scalar(@wgs_sample_list);
if ($wgs_sample_count > 0){
  open (OUT1, ">$vaf_file_matrix_wgs") || die "\n\nCould not open file: $vaf_file_matrix_wgs\n\n";
  open (OUT2, ">$mutation_status_file_matrix_wgs") || die "\n\nCould not open file: $mutation_status_file_matrix_wgs\n\n";
  my $wgs_sample_list_string = join("\t", @wgs_sample_list);
  my $wgs_header = "coord\t$wgs_sample_list_string";
  print OUT1 "$wgs_header\n";
  print OUT2 "$wgs_header\n";
  foreach my $coord (sort keys %rc){
    my @vafs;
    my @statuses;
    foreach my $sample (@wgs_sample_list){
      my $vaf = $rc{$coord}{$sample}{wgs_tumor_vaf};
      my $status = 0;
      if ($vaf > 0){
        $status = 1;
      }
      push(@vafs, $vaf);
      push(@statuses, $status);
    }
    my $wgs_val_string = join("\t", @vafs);
    my $wgs_status_string = join("\t", @statuses);
    print OUT1 "$coord\t$wgs_val_string\n";
    print OUT2 "$coord\t$wgs_status_string\n";
  }
  close(OUT1);
  close(OUT2);
}

#Exome VAFs
my @exome_sample_list = sort keys %exome_sample_list;
my $exome_sample_count = scalar(@exome_sample_list);
if ($exome_sample_count > 0){
  open (OUT1, ">$vaf_file_matrix_exome") || die "\n\nCould not open file: $vaf_file_matrix_exome\n\n";
  open (OUT2, ">$mutation_status_file_matrix_exome") || die "\n\nCould not open file: $mutation_status_file_matrix_exome\n\n";
  my $exome_sample_list_string = join("\t", @exome_sample_list);
  my $exome_header = "coord\t$exome_sample_list_string";
  print OUT1 "$exome_header\n";
  print OUT2 "$exome_header\n";
  foreach my $coord (sort keys %rc){
    my @vafs;
    my @statuses;
    foreach my $sample (@exome_sample_list){
      my $vaf = $rc{$coord}{$sample}{exome_tumor_vaf};
      my $status = 0;
      if ($vaf > 0){
        $status = 1;
      }
      push(@vafs, $vaf);
      push(@statuses, $status);
    }
    my $exome_val_string = join("\t", @vafs);
    my $exome_status_string = join("\t", @statuses);
    print OUT1 "$coord\t$exome_val_string\n";
    print OUT2 "$coord\t$exome_status_string\n";
  }
  close(OUT1);
  close(OUT2);
}

print BLUE, "\n\nWrote results to: $outdir\n\n", RESET;

exit();



