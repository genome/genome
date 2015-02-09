package Genome::Model::ClinSeq::Command::Converge::DgidbGenes;
#Written by Malachi Griffith

use strict;
use warnings;
use Genome;

class Genome::Model::ClinSeq::Command::Converge::DgidbGenes {
  is => ['Genome::Model::ClinSeq::Command::Converge::Base',
  'Genome::Model::ClinSeq::Util'],
  has_input => [
  outdir => {
    is => 'FilesystemPath',
    doc => 'Directory where output files will be written',
  },
  cancer_annotation_db => {
    is => 'Genome::Db::Tgi::CancerAnnotation',
    doc => 'cancer annotation db to be used.(not required if using a clinseq model as input)',
    example_values  => ['tgi/cancer-annotation/human/build37-20130401.1'],
  },
  ],
  doc => 'converge various types of Druggable genes results',
};

sub help_detail {
  return <<INFO
For a group of ClinSeq models, converge various types of Druggable genes results

Consider the following event types (individually, and then together):
SNVs, InDels, CNV amplifications, fusion genes, etc.
Create a summary at the gene level

Create a separate report for known druggable (e.g. DrugBank) and potentially druggable genes (e.g. kinases, etc.)
For potentially druggable summarize to Gene level as well as the Gene-Family level

Create a series of final reports that list genes/drugs and summarizes totals
 - genes with each kind of event
 - genes with any event
 - patients with at least 1 druggable gene (each type of event - then all event types together)
 - drugs indicated in at least 1 patient, 2 patients, etc.
 - total number of drugs
INFO
}

sub help_synopsis {
  return <<INFO;
  Example usage:
genome model clin-seq converge dgidb-genes --builds='id in ["4b7539bb10cc4b9c97577cf11f4c79a2","cdca0edf526c4fe193d3054627a5871b"]' --outdir=/tmp/snv_indel_report/

genome model clin-seq converge dgidb-genes --builds='model.model_groups.id=9d0fcdca2b5d4f4385b83d2f75addac4,is_last_complete=1' --outdir=/tmp/snv_indel_report/

genome model clin-seq converge dgidb-genes --builds='model_groups.id=9d0fcdca2b5d4f4385b83d2f75addac4,is_last_complete=1' --outdir=/tmp/snv_indel_report/

genome model clin-seq converge dgidb-genes --builds='model.id in ["279f50e35d2b479ea3c32486eafd4ad4","7143119a93984056ae3f32c88c9ac2a1"],is_last_complete=1' --outdir=/tmp/snv_indel_report/
INFO
}

sub execute {
  my $self = shift;
  my @clinseq_builds = $self->builds;
  my $cancer_annotation_db = $self->cancer_annotation_db;

  #Get files:
  my $files = $self->getFiles('-builds'=>\@clinseq_builds);
  my @event_types =  qw (snv indel cnv_gain rna_cufflinks_absolute rna_tophat_absolute);
  my $gene_symbol_lists_dir = File::Spec->catdir(
    $cancer_annotation_db->data_directory,
    "GeneSymbolLists/");
  my $gene_groups = 'Default';
  my $symbol_list_names = $self->importSymbolListNames('-gene_symbol_lists_dir'=>$gene_symbol_lists_dir, '-verbose'=>1);

  #Known druggable genes
  my $k_result = $self->parseKnownDruggableFiles('-files'=>$files, '-event_types'=>\@event_types);
  my $k_g = $k_result->{'genes'}; #Known druggable genes
  my $k_i = $k_result->{'interactions'}; #Known druggable gene interactions
  my $p_result = $self->parsePotentiallyDruggableFiles('-files'=>$files, '-event_types'=>\@event_types, '-symbol_list_names'=>$symbol_list_names, '-gene_groups'=>$gene_groups);
  my $p_gf = $p_result->{'gene_fams'}; #Potentially druggable gene families

  $self->generate_genefamily_patient_events($p_gf, \@event_types, $gene_groups,
    $symbol_list_names);
  $self->generate_interaction_patient_events($k_i, \@event_types);
  $self->generate_gene_individual_patients_events($files, $k_g, \@event_types);
  $self->generate_gene_patient_events($k_g, \@event_types);
  return 1;
}

sub generate_genefamily_patient_events {
  my $self = shift;
  my $p_gf = shift;
  my $event_types = shift;
  my $gene_groups = shift;
  my $symbol_list_names = shift;
  my $outdir = $self->outdir;

  #C.) Generate Gene Family -> patient events (either patient counts or hits) <- Gene Family Members (those that are actually affected)
  my @et_header = ();
  foreach my $et (@$event_types){
    push(@et_header, "$et"."_sample_count");
    push(@et_header, "$et"."_sample_list");
    push(@et_header, "$et"."_gene_count");
    push(@et_header, "$et"."_gene_list");
    push(@et_header, "$et"."_hit_count");
  }
  my $et_header_s = join("\t", @et_header);
  my $gene_family_table = $outdir . "/GeneFamilyTable.tsv";
  $self->status_message("Writing (gene family -> patient events <- gene family member) table to: $gene_family_table");
  open (OUT, ">$gene_family_table") || die "\nCould not open output file: $gene_family_table";
  print OUT "gene_group\tgene_group_size\tgrand_patient_count\tgrand_patient_list\tgrand_gene_count\tgrand_gene_list\tgrand_hit_count\t$et_header_s\n";
  my $sublists = $symbol_list_names->{sublists};
  my $target_groups = $sublists->{$gene_groups}->{groups};
  foreach my $group_name (sort {$target_groups->{$a}->{order} <=> $target_groups->{$b}->{order}} keys %{$target_groups}){
    my $group_gene_count = $p_gf->{$group_name}->{gene_count};
    unless (defined $p_gf->{$group_name}->{grand_list}) {
      next;
    }
    my %patient_list = %{$p_gf->{$group_name}->{grand_list}->{patient_list}};
    my @grand_patient_list = keys %patient_list;
    my @tmp = sort @grand_patient_list;
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
    foreach my $et (@$event_types){
      if (defined($p_gf->{$group_name}->{$et})){
        my %patients = %{$p_gf->{$group_name}->{$et}->{patient_list}};
        my @patient_list = keys %patients;
        my @tmp = sort @patient_list;;
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
}

sub generate_interaction_patient_events {
  my $self = shift;
  my $k_i = shift;
  my $event_types = shift;
  my $outdir = $self->outdir;
  my $known_interaction_table = $outdir . "/KnownInteractionTable.tsv";

  my @et_header;
  foreach my $et (@$event_types){
    push(@et_header, "$et"."_sample_count");
    push(@et_header, "$et"."_sample_list");
  }
  my $et_header_s = join("\t", @et_header);

  #B.) Generate Interaction -> patient events <- known drug table (only the one drug of the interaction)
  $self->status_message("Writing (interaction -> patient events <- known drug) table to: $known_interaction_table");
  open (OUT, ">$known_interaction_table") || die "\nCould not open output file: $known_interaction_table";
  print OUT "gene\tgrand_patient_count\tgrand_patient_list\tdrug_count\tdrug_list\t$et_header_s\n";
  foreach my $interaction (sort keys %{$k_i}){
    my $mapped_gene_name = $k_i->{$interaction}->{gene_name};
    my $grand_patient_list = $k_i->{$interaction}->{grand_list};
    my $grand_patient_count = keys %{$grand_patient_list};
    my @grand_patient_list = keys %{$grand_patient_list};
    my @tmp = sort @grand_patient_list;;
    my $grand_patient_list_s = join (",", @tmp);
    my $drug_list_s = $k_i->{$interaction}->{drug_name};

    my @et_values;
    foreach my $et (@$event_types){
      if (defined($k_i->{$interaction}->{$et})){
        my %patients = %{$k_i->{$interaction}->{$et}->{patient_list}};
        my @patient_list = keys %patients;
        my @tmp = sort @patient_list;;
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
}

sub generate_gene_individual_patients_events {
  my $self = shift;
  my $files = shift;
  my $k_g = shift;
  my $event_types = shift;
  my $outdir = $self->outdir;

  #A.) Generate Gene -> individual patients <- known drugs table (all drugs that target that gene)
  my @master_patient_list_u = keys %{$files};
  my @master_patient_list = sort @master_patient_list_u;
  my @pt_header;
  foreach my $pt (@master_patient_list){
    push(@pt_header, "$pt");
  }
  my $pt_header_s = join("\t", @pt_header);
  my $known_gene_drug_table2 = $outdir . "/KnownGeneDrugTable_IndividualPatients.tsv";
  $self->status_message("Writing (gene -> patient events <- known drugs) table to: $known_gene_drug_table2");
  open (OUT, ">$known_gene_drug_table2") || die "\nCould not open output file: $known_gene_drug_table2";
  print OUT "gene\tgrand_patient_count\tgrand_patient_list\tdrug_count\tdrug_list\tdrug_list_abr\tdrug_classes\t$pt_header_s\tsnv_aa_changes\tdb_source\tdb_filter\n";
  foreach my $gene (sort keys %{$k_g}){
    my $grand_patient_list = $k_g->{$gene}->{grand_list};
    my $grand_patient_count = keys %{$grand_patient_list};
    my @grand_patient_list = keys %{$grand_patient_list};
    my @tmp = sort @grand_patient_list;
    my $grand_patient_list_s = join (",", @tmp);
    my $drug_list = $k_g->{$gene}->{drug_list};
    my $drug_count = keys %{$drug_list};
    my @drug_list = keys %{$drug_list};
    my @sorted_drug_list = sort @drug_list;
    my $drug_list_s = '"' . join(",", @sorted_drug_list) . '"';
    my $drug_list_abr = $self->abbrevDrugList('-drug_list'=>\@sorted_drug_list);
    my $drug_class = $k_g->{$gene}->{drug_class};
    my @drug_class = keys %{$drug_class};
    my @sorted_drug_class = sort @drug_class;
    my $drug_class_s = join(",", @sorted_drug_class);
    my $aa_changes_s = "na";
    my $db_source = "na";
    my $filter = "na";
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
      foreach my $et (@$event_types){
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

    print OUT "$gene\t$grand_patient_count\t$grand_patient_list_s\t$drug_count\t$drug_list_s\t$drug_list_abr\t$drug_class_s\t$pt_values_s\t$aa_changes_s\t$db_source\t$filter\n";
  }
  close(OUT);
}

sub generate_gene_patient_events {
  my $self = shift;
  my $k_g = shift;
  my $event_types = shift;
  my $outdir = $self->outdir;
  #A.) Generate Gene -> patient events <- known drugs table
  # (all drugs that target that gene)
  #Generate a header for the event types output columns
  my @et_header;
  foreach my $et (@$event_types){
    push(@et_header, "$et"."_sample_count");
    push(@et_header, "$et"."_sample_list");
  }
  my $et_header_s = join("\t", @et_header);
  my $known_gene_drug_table1 = $outdir . "/KnownGeneDrugTable_PatientCounts.tsv";
  $self->status_message("Writing (gene -> patient events <- known drugs) table to: $known_gene_drug_table1");
  open (OUT, ">$known_gene_drug_table1") || die "\nCould not open output file: $known_gene_drug_table1";
  print OUT "gene\tgrand_patient_count\tgrand_patient_list\tdrug_count\tdrug_list\tdrug_list_abr\tdrug_classes\t$et_header_s\tsnv_aa_changes\tdb_source\tdb_filter\n";
  foreach my $gene (sort keys %{$k_g}){
    my $grand_patient_list = $k_g->{$gene}->{grand_list};
    my $grand_patient_count = keys %{$grand_patient_list};
    my @grand_patient_list = keys %{$grand_patient_list};
    my @tmp = sort @grand_patient_list;
    my $grand_patient_list_s = join (",", @tmp);
    my $drug_list = $k_g->{$gene}->{drug_list};
    my $drug_count = keys %{$drug_list};
    my @drug_list = keys %{$drug_list};
    my @sorted_drug_list = sort @drug_list;
    my $drug_list_s = join(",", @sorted_drug_list);
    my $drug_list_abr = $self->abbrevDrugList('-drug_list'=>\@sorted_drug_list);

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
    foreach my $et (@$event_types){
      if (defined($k_g->{$gene}->{$et})){
        my %patients = %{$k_g->{$gene}->{$et}->{patient_list}};
        my @patient_list = keys %patients;
        my @tmp = sort @patient_list;;
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
    print OUT "$gene\t$grand_patient_count\t$grand_patient_list_s\t$drug_count\t$drug_list_s\t$drug_list_abr\t$drug_class_s\t$et_values_s\t$aa_changes_s\n";
  }
  close(OUT);
}

#Abbreviate a long list of drugs (comma separated) into something a bit more readable
sub abbrevDrugList{
  my $self = shift;
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
  return('"' . $drug_list_abr . '"');
}

#Parse annotation files (containing drug gene family information)
sub parsePotentiallyDruggableFiles{
  my $self = shift;
  my %args = @_;
  my $files = $args{'-files'};
  my @event_types = @{$args{'-event_types'}};
  my $symbol_list_names = $args{'-symbol_list_names'};
  my $gene_groups = $args{'-gene_groups'};

  my $master_list = $symbol_list_names->{master_list};
  my $master_group_list = $symbol_list_names->{master_group_list};
  my $sublists = $symbol_list_names->{sublists};
  unless (defined($sublists->{$gene_groups})){
    die $self->error_message("Could not find a gene sublist with the name: $gene_groups");
  }
  my $target_groups = $sublists->{$gene_groups}->{groups};
  $self->status_message("Parsing files containing potentially druggable genes data");

  #Store all results organized by gene and separately by gene family (e.g. kinases, ion channels etc.)
  my %result;
  my %genes;
  my %gene_fams;

  #Note that the annotation files contain all qualifying variant events (the associated gene may or may not belong to a gene family of interest)
  foreach my $patient (keys %{$files}){
    $self->status_message("\n\t$patient");
    foreach my $event_type (@event_types){
      my $annot_file_path;
      foreach my $source (qw(wgs exome rnaseq)) {
        $annot_file_path = $files->{$patient}->{$event_type}->{$source}->{annot_file_path};
        if($annot_file_path) {
          last;
        }
      }
      unless($annot_file_path) {
        $self->warning_message("\nCould not find annotation file for $patient:$event_type");
        next;
      }
      $self->status_message("\n$patient\t$event_type\t\t$annot_file_path");
      open (IN, "$annot_file_path") || die "\nCould not open annotation file: $annot_file_path";
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
            unless ($value){
              next();
            }

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

1;

