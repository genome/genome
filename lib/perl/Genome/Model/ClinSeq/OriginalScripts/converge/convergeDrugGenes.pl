#!/usr/bin/env genome-perl
#Written by Malachi Griffith
#For a group of ClinSeq models, converge various types of Druggable genes results

#Consider the following event types (individually, and then together):
#SNVs, InDels, CNV amplifications, fusion genes, etc.
#Create a summary at the gene level

#Create a separate report for known druggable (e.g. DrugBank) and potentially druggable genes (e.g. kinases, etc.)
#For potentially druggable summarize to Gene level as well as the Gene-Family level

#Create a series of final reports that list genes/drugs and summarizes totals
# - genes with each kind of event
# - genes with any event
# - patients with at least 1 druggable gene (each type of event - then all event types together)
# - drugs indicated in at least 1 patient, 2 patients, etc.
# - total number of drugs

use strict;
use warnings;
use Getopt::Long;
use Term::ANSIColor qw(:constants);
use Data::Dumper;
use above 'Genome';
use Genome::Model::ClinSeq::Util qw(:all);
use Genome::Model::ClinSeq::OriginalScripts::Converge qw(:all);

my $reference_annotations_dir = '';
my $gene_groups = '';
my $build_ids = '';
my $model_ids = '';
my $model_group_id = '';
my $event_types_list = '';
my $dgidb_subdir_name = '';
my $filter_name = '';
my $outdir = '';
my $verbose = 0;

GetOptions ('reference_annotations_dir=s'=>\$reference_annotations_dir, 'gene_groups=s'=>\$gene_groups,
            'build_ids=s'=>\$build_ids, 'model_ids=s'=>\$model_ids, 'model_group_id=s'=>\$model_group_id,
            'event_types_list=s'=>\$event_types_list, 'dgidb_subdir_name=s'=>\$dgidb_subdir_name, 'filter_name=s'=>\$filter_name, 
            'outdir=s'=>\$outdir, 'verbose=i'=>\$verbose);

my $usage=<<INFO;
  Example usage: 

  convergeDrugGenes.pl  --model_group_id='31779'  --event_types_list='all'  --dgidb_subdir_name='drugbank'  --filter_name='antineo'  --outdir=/gscmnt/sata132/techd/mgriffit/luc/druggable_genes/  --reference_annotations_dir=/gscmnt/sata132/techd/mgriffit/reference_annotations/  --gene_groups='LUC17'  --verbose=1

  Specify the following to define which sets of gene groups will be summarized for potential druggability
  --reference_annotation_dirs  Path to reference annotation files
  --gene_groups                Name of gene group sublist (e.g. 'Default', 'LUC17').  See reference annotations dir for details

  Specify *one* of the following as input (each model/build should be a ClinSeq model)
  --build_ids                  Comma separated list of specific build IDs
  --model_ids                  Comma separated list of specific model IDs
  --model_group_id             A single genome model group ID

  Additional parameters
  --dgidb_subdir_name          The name of the drug-gene database source subdir (e.g. 'drugbank', 'santa_monica_lung')
  --filter_name                The name appended to each file indicating which filter was applied (e.g. 'default', 'antineo', etc.)
  --event_types_list           Specify a comma separated list of valid event types.  Use 'all' to include all types
                               Allowed event types: 'snv,indel,cnv_gain,rna_cufflinks_absolute,rna_tophat_absolute'

  --outdir                     Path of the a directory for output files
  --verbose                    More descriptive stdout messages

  Test Clinseq model groups:
  31779                        LUC17 project

INFO

unless ($reference_annotations_dir && $gene_groups && ($build_ids || $model_ids || $model_group_id) && $dgidb_subdir_name && $filter_name && $event_types_list && $outdir){
  print RED, "\n\nRequired parameter missing", RESET;
  print GREEN, "\n\n$usage", RESET;
  exit(1);
}

#Set output file names
$outdir = &checkDir('-dir'=>$outdir, '-clear'=>"no");
my $known_gene_drug_table1 = $outdir . "KnownGeneDrugTable_PatientCounts.tsv";
my $known_gene_drug_table2 = $outdir . "KnownGeneDrugTable_IndividualPatients.tsv";
my $known_interaction_table = $outdir . "KnownInteractionTable.tsv";
my $gene_family_table = $outdir . "GeneFamilyTable.tsv";

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

#Import gene symbol lists, groups of gene symbol lists, and the subset that is going to be processed for this analysis
my $gene_symbol_lists_dir = $reference_annotations_dir . "GeneSymbolLists/";
$gene_symbol_lists_dir = &checkDir('-dir'=>$gene_symbol_lists_dir, '-clear'=>"no");
my $symbol_list_names = &importSymbolListNames('-gene_symbol_lists_dir'=>$gene_symbol_lists_dir, '-verbose'=>$verbose);

#Get files:
# - drug-gene interaction files for each event type
# - annotated files for each event type containing potentially druggable results
my $files = &getFiles('-models_builds'=>$models_builds, '-filter_name'=>$filter_name, '-dgidb_subdir_name'=>$dgidb_subdir_name);

#Known druggable genes
my $k_result = &parseKnownDruggableFiles('-files'=>$files, '-event_types'=>\@event_types);
my $k_g = $k_result->{'genes'}; #Known druggable genes
my $k_i = $k_result->{'interactions'}; #Known druggable gene interactions
#print Dumper $k_g;

#Potentially druggable genes
my $p_result = &parsePotentiallyDruggableFiles('-files'=>$files, '-event_types'=>\@event_types, '-symbol_list_names'=>$symbol_list_names, '-gene_groups'=>$gene_groups);
my $p_gf = $p_result->{'gene_fams'}; #Potentially druggable gene families
#print Dumper $p_result; 


#Generate a header for the event types output columns
my @et_header;
foreach my $et (@event_types){
	push(@et_header, "$et"."_sample_count");
	push(@et_header, "$et"."_sample_list");
}
my $et_header_s = join("\t", @et_header);


#A.) Generate Gene -> patient events <- known drugs table (all drugs that target that gene)
print BLUE, "\n\nWriting (gene -> patient events <- known drugs) table to: $known_gene_drug_table1", RESET;
open (OUT, ">$known_gene_drug_table1") || die "\n\nCould not open output file: $known_gene_drug_table1\n\n";
print OUT "gene\tgrand_patient_count\tgrand_patient_list\tdrug_count\tdrug_list\tdrug_list_abr\tdrug_classes\t$et_header_s\tsnv_aa_changes\tdb_source\tdb_filter\n";
foreach my $gene (sort keys %{$k_g}){
	my $grand_patient_list = $k_g->{$gene}->{grand_list};
	my $grand_patient_count = keys %{$grand_patient_list};
	my @grand_patient_list = keys %{$grand_patient_list};
	my @tmp = sort { substr($a, &lengthOfAlpha($a)) <=> substr($b, &lengthOfAlpha($b)) } @grand_patient_list;
	my $grand_patient_list_s = join (",", @tmp);
	my $drug_list = $k_g->{$gene}->{drug_list};
	my $drug_count = keys %{$drug_list};
	my @drug_list = keys %{$drug_list};
	my @sorted_drug_list = sort @drug_list;
  my $drug_list_s = join(",", @sorted_drug_list);
  my $drug_list_abr = &abbrevDrugList('-drug_list'=>\@sorted_drug_list);

  my $drug_class = $k_g->{$gene}->{drug_class};
	my @drug_class = keys %{$drug_class};
	my @sorted_drug_class = sort @drug_class;
  my $drug_class_s = join(",", @sorted_drug_class);

  my $aa_changes_s = "na";
  if (defined($k_g->{$gene}->{aa_changes})){
    my %aa_changes = %{$k_g->{$gene}->{aa_changes}};
    my @aa_changes = keys %aa_changes;
    $aa_changes_s = join("|", @aa_changes);
  }
  
	my @et_values;
	foreach my $et (@event_types){
		if (defined($k_g->{$gene}->{$et})){
			my %patients = %{$k_g->{$gene}->{$et}->{patient_list}};
			my @patient_list = keys %patients;
			my @tmp = sort { substr($a, &lengthOfAlpha($a)) <=> substr($b, &lengthOfAlpha($b)) } @patient_list;;
      my $patient_list_s = join (",", @tmp);
			my $patient_count = scalar(@patient_list);
			push(@et_values, $patient_count);
			push(@et_values, $patient_list_s);
		}else{
			push(@et_values, 0);
			push(@et_values, "NA");
		}
	}
	my $et_values_s = join("\t", @et_values);
	print OUT "$gene\t$grand_patient_count\t$grand_patient_list_s\t$drug_count\t$drug_list_s\t$drug_list_abr\t$drug_class_s\t$et_values_s\t$aa_changes_s\t$dgidb_subdir_name\t$filter_name\n";
}
close(OUT);


#A.) Generate Gene -> individual patients <- known drugs table (all drugs that target that gene)
my @master_patient_list_u = keys %{$files};
my @master_patient_list = sort @master_patient_list_u;
my @pt_header;
foreach my $pt (@master_patient_list){
	push(@pt_header, "$pt");
}
my $pt_header_s = join("\t", @pt_header);

print BLUE, "\n\nWriting (gene -> patient events <- known drugs) table to: $known_gene_drug_table2", RESET;
open (OUT, ">$known_gene_drug_table2") || die "\n\nCould not open output file: $known_gene_drug_table2\n\n";
print OUT "gene\tgrand_patient_count\tgrand_patient_list\tdrug_count\tdrug_list\tdrug_list_abr\tdrug_classes\t$pt_header_s\tsnv_aa_changes\tdb_source\tdb_filter\n";
foreach my $gene (sort keys %{$k_g}){
	my $grand_patient_list = $k_g->{$gene}->{grand_list};
	my $grand_patient_count = keys %{$grand_patient_list};
	my @grand_patient_list = keys %{$grand_patient_list};
	my @tmp = sort { substr($a, &lengthOfAlpha($a)) <=> substr($b, &lengthOfAlpha($b)) } @grand_patient_list;
	my $grand_patient_list_s = join (",", @tmp);
	my $drug_list = $k_g->{$gene}->{drug_list};
	my $drug_count = keys %{$drug_list};
	my @drug_list = keys %{$drug_list};
	my @sorted_drug_list = sort @drug_list;
  my $drug_list_s = join(",", @sorted_drug_list);
  my $drug_list_abr = &abbrevDrugList('-drug_list'=>\@sorted_drug_list);
  my $drug_class = $k_g->{$gene}->{drug_class};
	my @drug_class = keys %{$drug_class};
	my @sorted_drug_class = sort @drug_class;
  my $drug_class_s = join(",", @sorted_drug_class);
  my $aa_changes_s = "na";
  if (defined($k_g->{$gene}->{aa_changes})){
    my %aa_changes = %{$k_g->{$gene}->{aa_changes}};
    my @aa_changes = keys %aa_changes;
    $aa_changes_s = join("|", @aa_changes);
  }

  #Assign each patient a status according to which event type is present
  my @pt_values;
  foreach my $pt (@master_patient_list){
    my $pt_value;
    my $single_value;
    my %pt_event_types;
	  foreach my $et (@event_types){
      if (defined($k_g->{$gene}->{$et})){
        my %et_patient_list = %{$k_g->{$gene}->{$et}->{patient_list}};
        if ($et_patient_list{$pt}){
          $pt_event_types{$et} = 1;
          $single_value = $et;
        }
      }
    }
    my $pt_et_count = keys %pt_event_types;
    if ($pt_et_count == 0){
      $pt_value = "NA";
    }elsif ($pt_et_count == 1){
      $pt_value = $single_value;
    }else{
      $pt_value = "multiple";
    }
    push (@pt_values, $pt_value);
  }
  my $pt_values_s = join("\t", @pt_values);

  print OUT "$gene\t$grand_patient_count\t$grand_patient_list_s\t$drug_count\t$drug_list_s\t$drug_list_abr\t$drug_class_s\t$pt_values_s\t$aa_changes_s\t$dgidb_subdir_name\t$filter_name\n";
}
close(OUT);



#B.) Generate Interaction -> patient events <- known drug table (only the one drug of the interaction)
print BLUE, "\n\nWriting (interaction -> patient events <- known drug) table to: $known_interaction_table", RESET;
open (OUT, ">$known_interaction_table") || die "\n\nCould not open output file: $known_interaction_table\n\n";
print OUT "gene\tgrand_patient_count\tgrand_patient_list\tdrug_count\tdrug_list\t$et_header_s\n";
foreach my $interaction (sort keys %{$k_i}){
  my $mapped_gene_name = $k_i->{$interaction}->{mapped_gene_name};
	my $grand_patient_list = $k_i->{$interaction}->{grand_list};
	my $grand_patient_count = keys %{$grand_patient_list};
	my @grand_patient_list = keys %{$grand_patient_list};
	my @tmp = sort { substr($a, &lengthOfAlpha($a)) <=> substr($b, &lengthOfAlpha($b)) } @grand_patient_list;;
	my $grand_patient_list_s = join (",", @tmp);
	my $drug_list_s = $k_i->{$interaction}->{drug_name};

	my @et_values;
	foreach my $et (@event_types){
		if (defined($k_i->{$interaction}->{$et})){
			my %patients = %{$k_i->{$interaction}->{$et}->{patient_list}};
			my @patient_list = keys %patients;
			my @tmp = sort { substr($a, &lengthOfAlpha($a)) <=> substr($b, &lengthOfAlpha($b)) } @patient_list;;
      my $patient_list_s = join (",", @tmp);
			my $patient_count = scalar(@patient_list);
			push(@et_values, $patient_count);
			push(@et_values, $patient_list_s);
		}else{
			push(@et_values, 0);
			push(@et_values, "NA");
		}
	}
	my $et_values_s = join("\t", @et_values);
	print OUT "$mapped_gene_name\t$grand_patient_count\t$grand_patient_list_s\t1\t$drug_list_s\t$et_values_s\n";
}
close(OUT);


#C.) Generate Gene Family -> patient events (either patient counts or hits) <- Gene Family Members (those that are actually affected)
@et_header = ();
foreach my $et (@event_types){
	push(@et_header, "$et"."_sample_count");
	push(@et_header, "$et"."_sample_list");
	push(@et_header, "$et"."_gene_count");
	push(@et_header, "$et"."_gene_list");
	push(@et_header, "$et"."_hit_count");
}
$et_header_s = join("\t", @et_header);

print BLUE, "\n\nWriting (gene family -> patient events <- gene family member) table to: $gene_family_table", RESET;
open (OUT, ">$gene_family_table") || die "\n\nCould not open output file: $gene_family_table\n\n";
print OUT "gene_group\tgene_group_size\tgrand_patient_count\tgrand_patient_list\tgrand_gene_count\tgrand_gene_list\tgrand_hit_count\t$et_header_s\n";


my $sublists = $symbol_list_names->{sublists};
my $target_groups = $sublists->{$gene_groups}->{groups};
foreach my $group_name (sort {$target_groups->{$a}->{order} <=> $target_groups->{$b}->{order}} keys %{$target_groups}){

  my $group_gene_count = $p_gf->{$group_name}->{gene_count};

  my %patient_list = %{$p_gf->{$group_name}->{grand_list}->{patient_list}};
  my @grand_patient_list = keys %patient_list;
	my @tmp = sort { substr($a, &lengthOfAlpha($a)) <=> substr($b, &lengthOfAlpha($b)) } @grand_patient_list;
  my $grand_patient_list_s = join (",", @tmp);
	my $grand_patient_count = scalar(@grand_patient_list);

  my %gene_list = %{$p_gf->{$group_name}->{grand_list}->{gene_list}};
  my @grand_gene_list = keys %gene_list;
  @tmp = sort @grand_gene_list;
  my $grand_gene_list_s = join (",", @tmp);
	my $grand_gene_count = scalar(@grand_gene_list);

  my $grand_hit_count = 0;
  foreach my $g (keys %gene_list){
    $grand_hit_count += $gene_list{$g};
  }

	my @et_values;
	foreach my $et (@event_types){
		if (defined($p_gf->{$group_name}->{$et})){
			my %patients = %{$p_gf->{$group_name}->{$et}->{patient_list}};
			my @patient_list = keys %patients;
			my @tmp = sort { substr($a, &lengthOfAlpha($a)) <=> substr($b, &lengthOfAlpha($b)) } @patient_list;;
      my $patient_list_s = join (",", @tmp);
			my $patient_count = scalar(@patient_list);
			push(@et_values, $patient_count);
			push(@et_values, $patient_list_s);

			my %genes = %{$p_gf->{$group_name}->{$et}->{gene_list}};
			my @gene_list = keys %genes;
      @tmp = sort @gene_list;
      my $gene_list_s = join (",", @tmp);
			my $gene_count = scalar(@gene_list);
			push(@et_values, $gene_count);
			push(@et_values, $gene_list_s);

      my $hit_count = 0;
      foreach my $g (keys %genes){
        $hit_count += $genes{$g};
      }
      push(@et_values, $hit_count);
		}else{
			push(@et_values, 0);
      push(@et_values, "NA");
			push(@et_values, 0);
      push(@et_values, "NA");
			push(@et_values, 0);
		}
	}
	my $et_values_s = join("\t", @et_values);
  print OUT "$group_name\t$group_gene_count\t$grand_patient_count\t$grand_patient_list_s\t$grand_gene_count\t$grand_gene_list_s\t$grand_hit_count\t$et_values_s\n";

}


close(OUT);


print "\n\n";

exit();




############################################################################################################################
#Abbreviate a long list of drugs (comma separated) into something a bit more readable                                      #
############################################################################################################################
sub abbrevDrugList{
  my %args = @_;
  my @drug_list = @{$args{'-drug_list'}};
  
  my $drug_list_abr = '';

  my $drug_list_count = scalar(@drug_list);

  if ($drug_list_count <= 3){
    $drug_list_abr = join(", ", @drug_list);
    $drug_list_abr .= ". [$drug_list_count]";
  }else{
    $drug_list_abr = "$drug_list[0], $drug_list[1], $drug_list[2], etc.";
    $drug_list_abr .= " [$drug_list_count]";
  }
  return($drug_list_abr);
}


############################################################################################################################
#Determine length of non-numeric portion of an alphanumeric string                                                         #
############################################################################################################################
sub lengthOfAlpha{
	my ($input) = @_;
	my ($alpha) = $input =~ /(\D{1,})/;
	my $length = length($alpha);
	return ($length);
}


############################################################################################################################
#Get input files to be parsed                                                                                              #
############################################################################################################################
sub getFiles{
  my %args = @_;
  my $models_builds = $args{'-models_builds'};
  my $filter_name = $args{'-filter_name'};
  my $dgidb_subdir_name = $args{'-dgidb_subdir_name'};
  my %files;

  if ($verbose){print BLUE, "\n\nGet annotation files and drug-gene interaction files from these builds", RESET;}
  my %mb = %{$models_builds->{cases}};
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
    my @snv_subdir_options = qw (wgs_exome wgs exome);
    my $snv_annot_file_name = "snvs.hq.tier1.v1.annotated.compact.tsv";
    my $snv_drug_file_name = "snvs.hq.tier1.v1.annotated.compact."."$filter_name".".tsv";
    foreach my $dir_name (@snv_subdir_options){
      my $annot_file_path = $topdir . "snv/$dir_name/$snv_annot_file_name";
      my $drug_file_path = $topdir . "snv/$dir_name/dgidb/$dgidb_subdir_name/$snv_drug_file_name";
      #If a file was already found, do nothing
      unless (defined($files{$final_name}{snv}{annot_file_path})){
        #If both files are present, store for later
        if (-e $annot_file_path && -e $drug_file_path){
          $files{$final_name}{snv}{annot_file_path} = $annot_file_path;
          $files{$final_name}{snv}{drug_file_path} = $drug_file_path;
        }
      }
    }
    #Make sure at least one pair of files was found
    unless (defined($files{$final_name}{snv}{annot_file_path}) && defined($files{$final_name}{snv}{drug_file_path})){
      print RED, "\n\nCould not find SNV drug-gene and annotation files for $final_name ($subject_name - $subject_common_name) in:\n\t$build_directory\n\n", RESET;
      exit(1);
    }

    #2.) Look for InDel files
    my @indel_subdir_options = qw (wgs_exome wgs exome);
    my $indel_annot_file_name = "indels.hq.tier1.v1.annotated.compact.tsv";
    my $indel_drug_file_name = "indels.hq.tier1.v1.annotated.compact."."$filter_name".".tsv";
    foreach my $dir_name (@indel_subdir_options){
      my $annot_file_path = $topdir . "indel/$dir_name/$indel_annot_file_name";
      my $drug_file_path = $topdir . "indel/$dir_name/dgidb/$dgidb_subdir_name/$indel_drug_file_name";
      #If a file was already found, do nothing
      unless (defined($files{$final_name}{indel}{annot_file_path})){
        #If both files are present, store for later
        if (-e $annot_file_path && -e $drug_file_path){
          $files{$final_name}{indel}{annot_file_path} = $annot_file_path;
          $files{$final_name}{indel}{drug_file_path} = $drug_file_path;
        }
      }
    }
    #Make sure at least one was found
    unless (defined($files{$final_name}{indel}{annot_file_path}) && defined($files{$final_name}{indel}{drug_file_path})){
      print RED, "\n\nCould not find INDEL drug-gene and annotation files for $final_name ($subject_name - $subject_common_name) in:\n\t$build_directory\n\n", RESET;
      exit(1);
    }

    #3.) Look for CNV gain files
    my $cnv_gain_annot_file_name = "cnv.AllGenes_Ensembl58.amp.tsv";
    my $cnv_gain_drug_file_name = "cnv.AllGenes_Ensembl58.amp."."$filter_name".".tsv";

    my $annot_file_path1 = $topdir . "cnv/$cnv_gain_annot_file_name";
    my $drug_file_path1 = $topdir . "cnv/dgidb/$dgidb_subdir_name/$cnv_gain_drug_file_name";
    my $annot_file_path2 = $topdir . "cnv/cnview/$cnv_gain_annot_file_name";
    my $drug_file_path2 = $topdir . "cnv/cnview/dgidb/$dgidb_subdir_name/$cnv_gain_drug_file_name";

    if (-e $annot_file_path1 && -e $drug_file_path1){
      $files{$final_name}{cnv_gain}{annot_file_path} = $annot_file_path1;
      $files{$final_name}{cnv_gain}{drug_file_path} = $drug_file_path1;
    }elsif(-e $annot_file_path2 && -e $drug_file_path2){
      $files{$final_name}{cnv_gain}{annot_file_path} = $annot_file_path2;
      $files{$final_name}{cnv_gain}{drug_file_path} = $drug_file_path2;      
    }else{
      print RED, "\n\nCould not find CNV drug-gene and annotation files for $final_name ($subject_name - $subject_common_name) in:\n\t$build_directory\n\n", RESET;
      exit(1);
    }

    #4.) Look for Cufflinks RNAseq outlier expression files 
    my $rna_cufflinks_annot_file_name = "isoforms.merged.fpkm.expsort.top1percent.tsv";
    my $rna_cufflinks_drug_file_name = "isoforms.merged.fpkm.expsort.top1percent."."$filter_name".".tsv";

    my $annot_file_path = $topdir . "rnaseq/tumor/cufflinks_absolute/isoforms_merged/$rna_cufflinks_annot_file_name";
    my $drug_file_path = $topdir . "rnaseq/tumor/cufflinks_absolute/isoforms_merged/dgidb/$dgidb_subdir_name/$rna_cufflinks_drug_file_name";
    if (-e $annot_file_path && -e $drug_file_path){
      $files{$final_name}{rna_cufflinks_absolute}{annot_file_path} = $annot_file_path;
      $files{$final_name}{rna_cufflinks_absolute}{drug_file_path} = $drug_file_path;
    }else{
      print RED, "\n\nCould not find Cufflinks Absolute drug-gene and annotation files for $final_name ($subject_name - $subject_common_name) in:\n\t$build_directory\n\t$annot_file_path\n\t$drug_file_path\n\n", RESET;
      exit(1);
    }

    #5.) Look for Tophat junction RNAseq outlier expression files
    my $rna_tophat_annot_file_name = "Ensembl.Junction.GeneExpression.top1percent.tsv";
    my $rna_tophat_drug_file_name = "Ensembl.Junction.GeneExpression.top1percent."."$filter_name".".tsv";

    $annot_file_path = $topdir . "rnaseq/tumor/tophat_junctions_absolute/$rna_tophat_annot_file_name";
    $drug_file_path = $topdir . "rnaseq/tumor/tophat_junctions_absolute/dgidb/$dgidb_subdir_name/$rna_tophat_drug_file_name";
    if (-e $annot_file_path && -e $drug_file_path){
      $files{$final_name}{rna_tophat_absolute}{annot_file_path} = $annot_file_path;
      $files{$final_name}{rna_tophat_absolute}{drug_file_path} = $drug_file_path;
    }else{
      print RED, "\n\nCould not find Tophat Absolute drug-gene and annotation files for $final_name ($subject_name - $subject_common_name) in:\n\t$build_directory\n\t$annot_file_path\n\t$drug_file_path\n\n", RESET;
      exit(1);
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

	#Note that the druggable files are already filtered down to only the variant affected genes with a drug interaction
  #To get a sense of the total number of events will have to wait until the annotation files are being proccessed 

  foreach my $patient (keys %{$files}){
    print BLUE, "\n\t$patient", RESET;
    foreach my $event_type (@event_types){
      my $drug_file_path = $files->{$patient}->{$event_type}->{drug_file_path};
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

        #If the event type is SNV.  Capture the AA changes
        if ($event_type eq "snv"){
          my $aa_changes = $line[$columns{'aa_changes'}{position}];
          if (defined($genes{$mapped_gene_name}{aa_changes})){
            my $aa_ref = $genes{$mapped_gene_name}{aa_changes};
            $aa_ref->{$aa_changes} = 1;
          }else{
            my %aa;
            $aa{$aa_changes} = 1;
            $genes{$mapped_gene_name}{aa_changes} = \%aa;
          }
        }

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
  $result{genes} = \%genes;
  $result{interactions} = \%interactions;

  return(\%result);
}




############################################################################################################################
#Parse annotation files (containing drug gene family information)                                                          #
############################################################################################################################
sub parsePotentiallyDruggableFiles{
  my %args = @_;
  my $files = $args{'-files'};
  my @event_types = @{$args{'-event_types'}};
  my $symbol_list_names = $args{'-symbol_list_names'};
  my $gene_groups = $args{'-gene_groups'};

  my $master_list = $symbol_list_names->{master_list};
  my $master_group_list = $symbol_list_names->{master_group_list};
  my $sublists = $symbol_list_names->{sublists};
  unless (defined($sublists->{$gene_groups})){
    print RED, "\n\nCould not find a gene sublist with the name: $gene_groups\n\n", RESET;
    exit();
  }
  my $target_groups = $sublists->{$gene_groups}->{groups};

  print BLUE, "\n\nParsing files containing potentially druggable genes data", RESET;

  #Store all results organized by gene and separately by gene family (e.g. kinases, ion channels etc.)
  my %result;
  my %genes;
  my %gene_fams;

	#Note that the annotation files contain all qualifying variant events (the associated gene may or may not belong to a gene family of interest)
  foreach my $patient (keys %{$files}){
    print BLUE, "\n\t$patient", RESET;
    foreach my $event_type (@event_types){
      my $annot_file_path = $files->{$patient}->{$event_type}->{annot_file_path};
      print BLUE, "\n\t\t$annot_file_path", RESET;

      open (IN, "$annot_file_path") || die "\n\nCould not open annotation file: $annot_file_path\n\n";
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

          #Check for columns of data needed for the specified gene families.  Throw warning every time a requested column is missing
          foreach my $group_name (sort {$target_groups->{$a}->{order} <=> $target_groups->{$b}->{order}} keys %{$target_groups}){
            my @list = @{$master_group_list->{$group_name}->{members}};
            foreach my $member (@list){
              unless($columns{$member}){
                print RED, "\n\nCould not find column for gene group: $group_name, member: $member\n\n", RESET;
                exit(1);
              }
            }
          }
          next();
        }
        my $mapped_gene_name = $line[$columns{'mapped_gene_name'}{position}];
 
        #Go through each gene group and aggregate patient by event type
        foreach my $group_name (sort {$target_groups->{$a}->{order} <=> $target_groups->{$b}->{order}} keys %{$target_groups}){
           my @list = @{$master_group_list->{$group_name}->{members}};
           my $group_gene_count = $master_group_list->{$group_name}->{gene_count};
           $gene_fams{$group_name}{gene_count} = $group_gene_count;
           foreach my $member (@list){

            #Skip undefined group member values? - meaningless if we choke on missing column above...
            unless ($columns{$member}){next();}

            #Get the value, for this patient, and this event type, and this gene group member (e.g. does this patient have an SNV in this gene that is a kinase)
            #If the value is 0, skip
            my $value = $line[$columns{$member}{position}];
            unless ($value){next();}
            
            if (defined($gene_fams{$group_name}{$event_type})){
              my $patients = $gene_fams{$group_name}{$event_type}{patient_list};
              $patients->{$patient}++;
              my $genes = $gene_fams{$group_name}{$event_type}{gene_list};
              $genes->{$mapped_gene_name}++;
            }else{
              my %patients;
              $patients{$patient} = 1;
              $gene_fams{$group_name}{$event_type}{patient_list} = \%patients;
              my %genes;
              $genes{$mapped_gene_name} = 1;
              $gene_fams{$group_name}{$event_type}{gene_list} = \%genes;
            }

            #Create or update the grand list of patients with ANY events hitting this gene
            if (defined($gene_fams{$group_name}{grand_list})){
              my $patients = $gene_fams{$group_name}{grand_list}{patient_list};
              $patients->{$patient}++;
              my $genes = $gene_fams{$group_name}{grand_list}{gene_list};
              $genes->{$mapped_gene_name}++;
            }else{
              my %patients;
              $patients{$patient} = 1;
              $gene_fams{$group_name}{grand_list}{patient_list} = \%patients;
              my %genes;
              $genes{$mapped_gene_name} = 1;
              $gene_fams{$group_name}{grand_list}{gene_list} = \%genes;
            }
          }
        }
      }
    }
  }

  $result{genes} = \%genes;
  $result{gene_fams} = \%gene_fams;
  return(\%result);
}



