#!/usr/bin/env genome-perl
#Written by Malachi Griffith

use strict;
use warnings;
use Getopt::Long;
use Term::ANSIColor qw(:constants);
use Data::Dumper;
use above 'Genome';
use Genome::Model::ClinSeq::Util qw(:all);

my $script_dir;
use Cwd 'abs_path';
BEGIN{
  if (abs_path($0) =~ /(.*\/).*\/.*\/.*\.pl/){
    $script_dir = $1;
  }
}

#Options / parameters
my $model_group_file = '';
my $working_dir = '';

GetOptions ('model_group_file=s'=>\$model_group_file, 'working_dir=s'=>\$working_dir);

my $usage=<<INFO;
  Example usage: 
  
  runDrugGeneInteractionComparisons.pl  --model_group_file=/gscmnt/sata132/techd/mgriffit/pcgp/DruggableGenes/SomaticModelGroupList.tsv  --working_dir=/gscmnt/sata132/techd/mgriffit/pcgp/DruggableGenes/results/

  Finds SNV files, calculates fixed gene names, adds info on model-groups, runs drug-gene-interaction code

  --model_group_file    File containing group names, case counts, and model-group IDs for somatic variation models
                        Format:  Group\tCaseCount\tModelGroupId

  --working_dir         Path to results directory

INFO

unless ($model_group_file && $working_dir){
  print RED, "\n\nRequired parameter missing", RESET;
  print GREEN, "\n\n$usage", RESET;
  exit(1);
}

#Check the working dir and create a subdirectory to store copies of the SNV files
$working_dir = &checkDir('-dir'=>$working_dir, '-clear'=>"no");
my $snv_subdir = $working_dir . "snv_files/";
mkdir ($snv_subdir);

#Set some output files
my $snvs_compact_outfile_genelevel = $working_dir . "all_cases_compact_snvs.tsv";
my $snvs_compact_outfile_genelevel_dgi_antineo = $working_dir . "all_cases_compact_snvs.antineoplastic.tsv";
my $snvs_compact_outfile_genelevel_dgi_inhibitor = $working_dir . "all_cases_compact_snvs.inhibitor.tsv";
my $snvs_compact_outfile_genelevel_dgi_kinase = $working_dir . "all_cases_compact_snvs.kinase.tsv";

#Get Entrez and Ensembl data for gene name mappings
my $entrez_ensembl_data = &loadEntrezEnsemblData();

#Directory of gene lists for various purposes
my $reference_annotations_dir = "/gscmnt/sata132/techd/mgriffit/reference_annotations/";
my $gene_symbol_lists_dir = $reference_annotations_dir . "GeneSymbolLists/";
$gene_symbol_lists_dir = &checkDir('-dir'=>$gene_symbol_lists_dir, '-clear'=>"no");
my @symbol_list_names1 = qw (Kinases KinasesGO CancerGeneCensus DrugBankAntineoplastic DrugBankInhibitors Druggable_RussLampel TfcatTransFactors FactorBookTransFactors TranscriptionFactorBinding_GO0008134 TranscriptionFactorComplex_GO0005667 CellSurface_GO0009986 DnaRepair_GO0006281 DrugMetabolism_GO0017144 TransporterActivity_GO0005215 ExternalSideOfPlasmaMembrane_GO0009897 GpcrActivity_GO0045028 GrowthFactorActivity_GO0008083 HistoneModification_GO0016570 HormoneActivity_GO0005179 IonChannelActivity_GO0005216 LipidKinaseActivity_GO0001727 NuclearHormoneReceptor_GO0004879 PeptidaseInhibitorActivity_GO0030414 PhospholipaseActivity_GO0004620 PhospoproteinPhosphataseActivity_GO0004721 ProteinSerineThreonineKinaseActivity_GO0004674 ProteinTyrosineKinaseActivity_GO0004713 RegulationOfCellCycle_GO0051726 ResponseToDrug_GO0042493);
my $gene_symbol_lists = &importGeneSymbolLists('-gene_symbol_lists_dir'=>$gene_symbol_lists_dir, '-symbol_list_names'=>\@symbol_list_names1, '-entrez_ensembl_data'=>$entrez_ensembl_data, '-verbose'=>0);


#Get list of model group IDs and associate group names
print BLUE, "\n\nImporting list of model groups from: $model_group_file", RESET;
my %mg_list;
open (MGLIST, "$model_group_file") || die "\n\nCould not open model-group file: $model_group_file\n\n";
my $header = 1;
while(<MGLIST>){
  chomp($_);
  my @line = split("\t", $_);
  if ($header){
    $header = 0;
    next();
  }
  my $group_name = $line[0];
  my $case_count = $line[1];
  my $model_group_id = $line[2];
  $mg_list{$group_name}{id} = $model_group_id;
  $mg_list{$group_name}{case_count} = $case_count;
}
close(MGLIST);
my $group_count = keys %mg_list;

#Summarize at the Gene and SNV level
my %master_snvs;
my %master_genes;
my %patient_list;
my %sample_list;

print BLUE, "\n\nGetting build directories and associated model/build info from the genome API for $group_count Model-Groups", RESET;
my $c = 0;
foreach my $group_name (sort keys %mg_list){
  my $mg_id = $mg_list{$group_name}{id};
  my $case_count = $mg_list{$group_name}{case_count};

  #Get model group object
  #genome/lib/perl/Genome/ModelGroup.pm
  my $mg = Genome::ModelGroup->get("id"=>$mg_id);

  #Get display name of the model group
  my $display_name = $mg->__display_name__;

  #Get subjects (e.g. samples) associated with each model of the model-group
  my @subjects = $mg->subjects;

  #Get the members of the model-group, i.e. the models
  my @models = $mg->models;

  #Cycle through the models and get their builds
  foreach my $m (@models){
    $c++;
    my $model_name = $m->name;
    my $model_id = $m->genome_model_id;
    my $last_complete_build_id = $m->_last_complete_build_id || "NONE_FINISHED";
    my $individual_id = $m->subject->patient->id;
    my $subject_name = $m->subject->patient->name;
    my $common_name = $m->subject->patient->common_name;
    my $sample_id = $m->subject->id;
    my $sample = Genome::Sample->get($sample_id);

    my $sample_name = $sample->name || "UNDEF";
    my $extraction_type = $sample->extraction_type || "UNDEF";
    my $sample_common_name = $sample->common_name || "UNDEF";
    my $tissue_desc = $sample->tissue_desc || "UNDEF";
    my $cell_type = $sample->cell_type || "UNDEF";
    my $effects_directory = "NA";
    if ($m->last_succeeded_build){
       my $b = $m->last_succeeded_build;
       $effects_directory = $b->data_directory . "/effects/";
    }

    #Make sure the SNV files exist
    my $snv_anno_top_file_name = "snvs.hq.tier1.v1.annotated.top";
    my $snv_anno_top_file = $effects_directory . $snv_anno_top_file_name;
    my $snv_anno_file_name = "snvs.hq.tier1.v1.annotated";
    my $snv_anno_file = $effects_directory . $snv_anno_file_name;

    #Some basic checks
    unless ($m->last_succeeded_build){
      print YELLOW, "\n\tCould not find a successful build for model ID: $model_id", RESET;
      next();
    }
    unless (-e $snv_anno_file){
      print YELLOW, "\n\tCould not find an SNV annotated file ($snv_anno_file)", RESET;
      next();
    }
    print BLUE, "\n\t$c\t$common_name\t$group_name\t$sample_common_name\t$sample_name\t$effects_directory", RESET;
    $patient_list{$common_name}=1;
    $sample_list{$sample_name}=1;

    #Make a copy of the SNV anno file to a subdir in the working dir
    my $new_file = $snv_subdir . "$sample_name"."_"."$snv_anno_file_name";
    my $cp_cmd = "cp $snv_anno_file $new_file";
    Genome::Sys->shellcmd(cmd => $cp_cmd);

    #Grab the SNV records - transcript level and collapse multiple transcript records to one record per SNV position
    open (SNV_IN, "$snv_anno_file") || die "\n\nCould not open SNV input file: $snv_anno_file\n\n";
    my %snvs;
    while(<SNV_IN>){
      chomp($_);
      my @line = split("\t", $_);
      my $chr = $line[0];
      my $coord = "$line[0]:$line[1]-$line[2]";
      my $gene_name = $line[6];
      my $mutation_type = $line[13];
      my $aa_change = $line[15];

      #Skip MT mutation
      if ($chr =~ /MT/i){
        next();
      }
      #Only allow SNVs of the type: missense, nonsense, splice_site, rna, splice_region
      my $snv_filter = "missense|nonsense|splice_site|splice_region|rna";
      unless ($mutation_type =~ /$snv_filter/){
        next();
      }
      my $mapped_gene_name = &fixGeneName('-gene'=>$gene_name, '-entrez_ensembl_data'=>$entrez_ensembl_data, '-verbose'=>0);
      if ($snvs{$coord}){
        my $mut_types = $snvs{$coord}{mut_types};
        $mut_types->{$mutation_type} = 1;
        my $aa_changes = $snvs{$coord}{aa_changes};
        $aa_changes->{$aa_change} = 1;
      }else{
        $snvs{$coord}{gene_name} = $gene_name;
        $snvs{$coord}{mapped_gene_name} = $mapped_gene_name;
        my %mut_types;
        $mut_types{$mutation_type} = 1;
        $snvs{$coord}{mut_types} = \%mut_types;
        my %aa_changes;
        $aa_changes{$aa_change} = 1;
        $snvs{$coord}{aa_changes} = \%aa_changes;
      }
      #print "\n\tDEBUG: $coord\t$gene_name\t$mapped_gene_name\t$mutation_type\t$aa_change";
    }
    close(SNV_IN);

    #Now go through the SNVs and store at the gene level - using the original gene name as a key
    foreach my $coord (sort keys %snvs){
      my $gene_name = $snvs{$coord}{gene_name};
      my $mapped_gene_name = $snvs{$coord}{mapped_gene_name};
      my $aa_changes = $snvs{$coord}{aa_changes};
      my $mut_types = $snvs{$coord}{mut_types};

      if ($master_genes{$gene_name}){

        #AA changes merged to the gene level (retaining transcript specific aa changes)
        my $gene_aa_changes = $master_genes{$gene_name}{gene_aa_changes};
        foreach my $aa_change (keys %{$aa_changes}){
          if ($gene_aa_changes->{$aa_change}){
            $gene_aa_changes->{$aa_change}->{count}++;
          }else{
            $gene_aa_changes->{$aa_change}->{count}=1;
          }
        }
        #Mutation types merged to the gene level (retaining transcript specific mutation_types)
        my $gene_mut_types = $master_genes{$gene_name}{gene_mut_types};
        foreach my $mut_type (keys %{$mut_types}){
          if ($gene_mut_types->{$mut_type}){
            $gene_mut_types->{$mut_type}++;
          }else{
            $gene_mut_types->{$mut_type}=1;
          }
        }
        
        #Master list of distinct mutation positions for this gene
        my $positions = $master_genes{$gene_name}{mutant_positions};
        if ($positions->{$coord}){
          $positions->{$coord}++;          
        }else{
          $positions->{$coord} = 1;
        }

        #Master list of common names (patient cases) having a mutation in this gene
        my $common_names = $master_genes{$gene_name}{common_names};
        $common_names->{$common_name} = 1;

        #Master list of cancer group names (cancer types) having a mutation in this gene
        my $group_names = $master_genes{$gene_name}{group_names};
        $group_names->{$group_name} = 1;

      }else{
        $master_genes{$gene_name}{mapped_gene_name} = $mapped_gene_name;

        #AA changes merged to the gene level (retaining transcript specific aa changes)
        my %gene_aa_changes;
        foreach my $aa_change (keys %{$aa_changes}){
          $gene_aa_changes{$aa_change}{count} = 1;
        }
        $master_genes{$gene_name}{gene_aa_changes} = \%gene_aa_changes;

        #Mutation types merged to the gene level (retaining transcript specific mutation types)
        my %gene_mut_types;
        foreach my $mut_type (keys %{$mut_types}){
          $gene_mut_types{$mut_type} = 1;
        }
        $master_genes{$gene_name}{gene_mut_types} = \%gene_mut_types;

        #Master list of distinct mutation positions for this gene
        my %positions;
        $positions{$coord} = 1;
        $master_genes{$gene_name}{mutant_positions} = \%positions;

        #Master list of common names (patient cases) having a mutation in this gene
        my %common_names;
        $common_names{$common_name} = 1;
        $master_genes{$gene_name}{common_names} = \%common_names;

        #Master list of cancer group names (cancer types) having a mutation in this gene
        my %group_names;
        $group_names{$group_name} = 1;
        $master_genes{$gene_name}{group_names} = \%group_names;
      }
    }
  }
}

my $patient_count = keys %patient_list;
my $sample_count = keys %sample_list;
print BLUE, "\n\nImported data for: $patient_count patients and $sample_count samples", RESET;

#Print out the gene level summary file
open (GENE_OUT, ">$snvs_compact_outfile_genelevel") || die "\n\nCould not open output file: $snvs_compact_outfile_genelevel\n\n";
print GENE_OUT "gene_name\tmapped_gene_name\thotspot_ratio\tmax_position_recurrence\tmutant_cases\tmutant_case_count\tgroup_names_string\tgroup_names_count\tmutant_positions\tmutant_positions_count\tmutation_types\tamino_acid_changes\n";
foreach my $gene_name (sort keys %master_genes){
  my $mapped_gene_name = $master_genes{$gene_name}{mapped_gene_name};
  my $gene_aa_changes = $master_genes{$gene_name}{gene_aa_changes};

  my $aa_change_string = '';
  foreach my $aa_change (sort {$gene_aa_changes->{$b}->{count} <=> $gene_aa_changes->{$a}->{count}} keys %{$gene_aa_changes}){
    $aa_change_string .= "$aa_change [$gene_aa_changes->{$aa_change}->{count}], ";
  }

  my $gene_mut_types = $master_genes{$gene_name}{gene_mut_types};
  my @gene_mut_types = keys %{$gene_mut_types};
  my @gene_mut_types_sort = sort @gene_mut_types;
  my $gene_mut_types_string = join(",", @gene_mut_types_sort);

  my $positions = $master_genes{$gene_name}{mutant_positions};
  my $positions_count = keys %{$positions};
  my @pos = keys %{$positions};
  my @pos_sort = sort @pos;
  my $pos_string = join (",", @pos_sort);

  my $max_position_recurrence = 0;
  foreach my $pos (keys %{$positions}){
    my $count = $positions->{$pos};
    if ($count > $max_position_recurrence){
      $max_position_recurrence = $count;
    }
  }
  my $common_names = $master_genes{$gene_name}{common_names};
  my $common_names_count = keys %{$common_names};
  my @common_names = keys %{$common_names};
  my @common_names_sort = sort @common_names;
  my $common_names_string = join (",", @common_names_sort);
  my $hotspot_ratio = $common_names_count / $positions_count;

  my $group_names = $master_genes{$gene_name}{group_names};
  my $group_names_count = keys %{$group_names};
  my @group_names = keys %{$group_names};
  my @group_names_sort = sort @group_names;
  my $group_names_string = join (",", @group_names_sort);

  print GENE_OUT "$gene_name\t$mapped_gene_name\t$hotspot_ratio\t$max_position_recurrence\t$common_names_string\t$common_names_count\t$group_names_string\t$group_names_count\t$pos_string\t$positions_count\t$gene_mut_types_string\t$aa_change_string\n";
}
close(GENE_OUT);

#print Dumper %master_genes;

#Run the drug-gene interaction script on the output file just created
my $dgi_cmd = "$script_dir"."summary/identifyDruggableGenes.pl  --candidates_file=$snvs_compact_outfile_genelevel --name_col_1=2  --interactions_file=/gscmnt/sata132/techd/mgriffit/DruggableGenes/KnownDruggable/DrugBank/query_files/DrugBank_WashU_INTERACTIONS.filtered.4.tsv  --name_col_2=12 > $snvs_compact_outfile_genelevel_dgi_antineo";
print YELLOW, "\n\n$dgi_cmd", RESET;
Genome::Sys->shellcmd(cmd => $dgi_cmd);

$dgi_cmd = "$script_dir"."summary/identifyDruggableGenes.pl  --candidates_file=$snvs_compact_outfile_genelevel --name_col_1=2  --interactions_file=/gscmnt/sata132/techd/mgriffit/DruggableGenes/KnownDruggable/DrugBank/query_files/DrugBank_WashU_INTERACTIONS.filtered.5.tsv  --name_col_2=12 > $snvs_compact_outfile_genelevel_dgi_inhibitor";
print YELLOW, "\n\n$dgi_cmd", RESET;
Genome::Sys->shellcmd(cmd => $dgi_cmd);

$dgi_cmd = "$script_dir"."summary/identifyDruggableGenes.pl  --candidates_file=$snvs_compact_outfile_genelevel --name_col_1=2  --interactions_file=/gscmnt/sata132/techd/mgriffit/DruggableGenes/KnownDruggable/DrugBank/query_files/DrugBank_WashU_INTERACTIONS.filtered.6.tsv  --name_col_2=12 > $snvs_compact_outfile_genelevel_dgi_kinase";
print YELLOW, "\n\n$dgi_cmd", RESET;
Genome::Sys->shellcmd(cmd => $dgi_cmd);

#Parse the antineoplastic gene hits to allow filtering of the file below?
open (DRUGGABLE, "$snvs_compact_outfile_genelevel_dgi_antineo") || die "\n\nCould not open file: $snvs_compact_outfile_genelevel_dgi_antineo\n\n";
my %druggable_genes;
$header = 1;
while (<DRUGGABLE>){
  chomp($_);
  my @line = split("\t", $_);
  if ($header){
    $header = 0;
    next();
  }
  $druggable_genes{$line[1]} = 1;
}
close(DRUGGABLE);

#Perform potentially druggable genes analysis by intersecting with various gene categories
my $path = $snvs_compact_outfile_genelevel;
print BLUE, "\n\tProcessing: $path", RESET;
my $new_path1 = $path.".gene.families.tsv";
my $new_path2 = $path.".gene.families.antineo.filtered.tsv";
open (INDATA, "$path") || die "\n\nCould not open input datafile: $path\n\n";
my %data;
my %cols;
$header = 1;
my $header_line = '';
my $l = 0;
while(<INDATA>){
  $l++;
  chomp($_);
  my $record = $_;
  my @line = split("\t", $_);       
  if ($header == 1){
    my $c = 0;
    $header_line = $_;
    foreach my $colname (@line){
      $cols{$colname}{position} = $c;
      $c++;
    }
    $header = 0;
    unless ($cols{'mapped_gene_name'}){
      print RED, "\n\nFile has no 'mapped_gene_name' column: $path\n\n", RESET;
      exit (1);
    }
    next();
  }
  $data{$l}{record} = $record;
  $data{$l}{mapped_gene_name} = $line[$cols{'mapped_gene_name'}{position}];
  $data{$l}{gene_name} = $line[$cols{'gene_name'}{position}];

}
close(INDATA);

#Figure out the gene matches to the gene symbol lists
#Test each gene name in this column against those in the list and add a column with the match status (i.e. is is a kinase, cancer gene, etc.)
foreach my $l (keys %data){
  my $mapped_gene_name = $data{$l}{mapped_gene_name};
  foreach my $gene_symbol_type (keys %{$gene_symbol_lists}){
    my $gene_symbols = $gene_symbol_lists->{$gene_symbol_type}->{symbols};
    if ($gene_symbols->{$mapped_gene_name}){
      $data{$l}{$gene_symbol_type} = 1;
    }else{
      $data{$l}{$gene_symbol_type} = 0;
    }
  }
}

#Print out a new file contain the extra columns
open (OUTDATA1, ">$new_path1") || die "\n\nCould not open output datafile: $new_path1\n\n";
open (OUTDATA2, ">$new_path2") || die "\n\nCould not open output datafile: $new_path2\n\n";
my @gene_symbol_list_names = sort {$gene_symbol_lists->{$a}->{order} <=> $gene_symbol_lists->{$b}->{order}} keys %{$gene_symbol_lists};
my $gene_symbol_list_name_string = join("\t", @gene_symbol_list_names);
print OUTDATA1 "$header_line\t$gene_symbol_list_name_string\n";
print OUTDATA2 "$header_line\t$gene_symbol_list_name_string\n";
foreach my $l (sort {$a <=> $b} keys %data){
  my $gene_name = $data{$l}{gene_name};
  my @tmp;
  foreach my $gene_symbol_list_name (@gene_symbol_list_names){
    push (@tmp, $data{$l}{$gene_symbol_list_name});
  }
  my $new_cols_string = join("\t", @tmp);
  print OUTDATA1 "$data{$l}{record}\t$new_cols_string\n";

  unless ($druggable_genes{$gene_name}){
    print OUTDATA2 "$data{$l}{record}\t$new_cols_string\n";
  }

}
close(OUTDATA1);
close(OUTDATA2);

print "\n\n";

exit();


