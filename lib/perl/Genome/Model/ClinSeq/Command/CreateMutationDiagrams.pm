package Genome::Model::ClinSeq::Command::CreateMutationDiagrams;

#Written by Malachi Griffith

use strict;
use warnings;
use Genome;
use Data::Dumper;
use Genome::Model::ClinSeq::Util qw(:all);

class Genome::Model::ClinSeq::Command::CreateMutationDiagrams {
    is => 'Command::V2',
    has_input => [
        cancer_annotation_db => {
            is => 'Genome::Db::Tgi::CancerAnnotation',
            doc => 'cancer-related gene lists, and associated data',
        },
        cosmic_annotation_db => {
            is => 'Genome::Db::Cosmic',
            doc => 'data from COSMIC (catalog of somatic mutations in cancer',
        },
        builds => { 
              is => 'Genome::Model::Build::SomaticVariation',
              is_many => 1,
              shell_args_position => 1,
              require_user_verify => 0,
              doc => 'somatic variation build(s) to create mutation diagrams from',
        },
        outdir => { 
              is => 'FilesystemPath',
              doc => 'Directory where output files will be written', 
        },
        collapse_variants => {
              is => 'Boolean',
              is_optional => 1,
              doc => 'Remove identical variants from different somatic-variation builds (do this for WGS + Exome builds from a single case)',
        },
        effect_type_filter => {
              is => 'Text',
              is_optional => 1,
              default => 'silent|3_prime_flanking_region|5_prime_flanking_region|intronic|3_prime_untranslated_region|5_prime_untranslated_region',
              doc => 'Translation effect types to filter out of results before generating mutation diagrams',
        },
        max_snvs_per_file => {
              is => 'Number',
              is_optional => 1,
              doc => 'Use this parameter to exclude SNV variant files that have a large number of variants (e.g. hyper-mutant cases)',
        },
        max_indels_per_file => {
              is => 'Number',
              is_optional => 1,
              doc => 'Use this parameter to exclude INDEL variant files that have a large number of variants (e.g. hyper-mutant cases)',
        },
        max_cosmic_sites => {
              is => 'Number',
              is_optional => 1,
              default => 25,
              doc => 'The max number of cosmic amino acid mutation positions to show for a single transcript',
        },
        max_count_per_cosmic_site => {
              is => 'Number',
              is_optional => 1,
              default => 10,
              doc => 'The max number of observations to display for a single COSMIC mutation',
        },
        max_transcripts => {
              is => 'Number',
              is_optional => 1,
              doc => 'For debugging purposes, limit to this number of transcripts',
        },
    ],
    doc => 'summarize the SVs of somatic variation build',
};

sub help_synopsis {
    return <<EOS

genome model clin-seq create-mutation-diagrams --outdir=/tmp/  129973671

genome model clin-seq create-mutation-diagrams --outdir=/tmp/  id=129973671

genome model clin-seq create-mutation-diagrams --outdir=/tmp/  model.id=2888915570

genome model clin-seq create-mutation-diagrams --outdir=/tmp/  "model.name='H_LF-09-213F-1221853.somatic_variation-2'"

genome model clin-seq create-mutation-diagrams --outdir=/tmp/  'id in [129973671,129708625]'

genome model clin-seq create-mutation-diagrams --outdir=/tmp/create_mutation_diagram/  --collapse-variants  --max-snvs-per-file=250  --max-indels-per-file=250  'id in [129973671,129708625]'

EOS
}

sub help_detail {
    return <<EOS
Create mutation diagrams ('loli-plots') for one or more somatic-variation builds.  Only annotated Tier1 variants will be used

A similar plot will be created for Cosmic mutations based on the annotation build of the somatic-variation builds specified

In clin-seq it would be typical to specify and *collapse* two somatic-variation builds, one from WGS and one Exome data

(put more content here)
EOS
}

sub __errors__ {
  my $self = shift;
  my @errors = $self->SUPER::__errors__(@_);

  unless (-e $self->outdir && -d $self->outdir) {
      push @errors, UR::Object::Tag->create(
	                                          type => 'error',
	                                          properties => ['outdir'],
	                                          desc => "Outdir: " . $self->outdir . " not found or not a directory",
                                          );
  }
  return @errors;
}

sub execute {
  my $self = shift;
  my @builds = $self->builds;
  my $outdir = $self->outdir;

  unless ($outdir =~ /\/$/){
    $outdir .= "/";
  }

  #Get the annotation reference name for the set of somatic-variation builds supplied
  my $annotation_reference_build = $self->resolve_annotation_build;

  #Get COSMIC mutation annotation for the annotation build used for this somatic variation model and cosmic version specified by the user
  #TODO: Create official versions of these data on allocated disk
  my $cancer_annotation_db = $self->cancer_annotation_db;
  my $cosmic_annotation_db = $self->cosmic_annotation_db;

  my $clinseq_annotations_dir = $cancer_annotation_db->data_directory . '/'; 
  my $ab_name = $annotation_reference_build->name;
  my $ab_name_f = $ab_name;
  $ab_name_f =~ s/\//\./g;


  my $cosmic_annotation_file = $cosmic_annotation_db->data_directory . "/clinseq/cosmic." . $ab_name_f . ".anno.filt";

  #Define final output files
  my $complete_somatic_variants_file = $outdir . "variants.hq.tier1.v1.annotated";
  my $filtered_somatic_variants_file = $outdir . "variants.hq.tier1.v1.annotated.filtered";
  my $filtered_cosmic_variants_file = $outdir . "cosmic.annotated.filtered";

  #Import Tier1 variants from the somatic-variation builds - key on transcript IDs
  my $somatic_variants = $self->import_somatic_variants('-complete_variants_file'=>$complete_somatic_variants_file, '-filtered_variants_file'=>$filtered_somatic_variants_file);

  #Import Cosmic SNVs from the Cosmic annotation file but only for those transcripts with imported Tier1 SNVs
  my $cosmic_variants = $self->import_cosmic_variants('-somatic_variants'=>$somatic_variants, '-cosmic_annotation_file'=>$cosmic_annotation_file, '-filtered_variants_file'=>$filtered_cosmic_variants_file);

  #Create mutation diagrams for every transcript (one for the somatic-variation data and one for the cosmic data)
  $self->draw_mutation_diagrams('-somatic_variants_file'=>$filtered_somatic_variants_file, '-cosmic_variants_file'=>$filtered_cosmic_variants_file, '-annotation_build_name'=>$ab_name);

  $self->status_message("\n\nCOMPLETE\n\n");
  return 1;
}


#####################################################################################################################
#Determine the reference annotation build from the input somatic variation builds                                   #
#####################################################################################################################
sub resolve_annotation_build{
  my $self = shift;
  my @builds = $self->builds;

  #Get the annotation version for the series of somatic variation builds specified
  my %ab_names;
  my %tmp;
  my $final_annotation_build;
  for my $somatic_build (@builds) {
    my $build_id = $somatic_build->id;
    my $model = $somatic_build->model;
    my $model_name = $model->name;
    my $annotation_build = $model->annotation_build;
    my $ab_name = $annotation_build->name;
    next unless $ab_name;
    $tmp{$ab_name}=1;
    $ab_names{$build_id}{annotation_build_name}=$ab_name;
    $ab_names{$build_id}{somatic_model_name}=$model_name;
    $final_annotation_build = $annotation_build;
  }
  if (keys %tmp > 1){
    $self->error_message("Found non-matching annotation build names for the list of somatic-variation builds");
    print Dumper %ab_names;
    exit 1;
  }
  unless ($final_annotation_build){
    $self->error_message("Unable to resolve annotation build for the list of somatic-variation builds");
    exit 1;
  }
  return $final_annotation_build;
}


#####################################################################################################################
#Import variants from the somatic-variation effects dirs                                                            #
#####################################################################################################################
sub import_somatic_variants{
  my $self = shift;
  my %args = @_;
  my $complete_variants_file = $args{'-complete_variants_file'};
  my $filtered_variants_file = $args{'-filtered_variants_file'};
  my @builds = $self->builds;

  my %var;

  #Hardcode location of files.  TODO: Add option to supply your own annotation files in case you want to do something custom...
  my @file_list;
  for my $somatic_build (@builds) {
    my $build_dir = $somatic_build->data_directory;
    my $tier1_snv_file = $build_dir . "/effects/snvs.hq.tier1.v1.annotated";
    my $tier1_indel_file = $build_dir . "/effects/indels.hq.tier1.v1.annotated";
    
    #Add the SNV file to the list to be processed - if an optional max-allowed SNVs parameter was supplied, determine the number of unique SNVs in this file and skip if neccessary
    unless (-e $tier1_snv_file){
      $self->error_message("\n\nCould not find tier1 SNV annotated file in build dir: $build_dir\n\n");
      exit 1;
    }
    if ($self->max_snvs_per_file){
      my $var_count = $self->unique_variant_count('-variant_file'=>$tier1_snv_file);
      if ($var_count > $self->max_snvs_per_file){
        $self->status_message("Too many variants ... skipping file: $tier1_snv_file");
      }else{
        push(@file_list, $tier1_snv_file);
      }
    }else{
      push(@file_list, $tier1_snv_file);
    }

    #Add the INDEL file to the list to be processed - if an optional max-allowed SNVs parameter was supplied, determine the number of unique SNVs in this file and skip if neccessary
    unless (-e $tier1_indel_file){
      $self->error_message("\n\nCould not find tier1 INDEL annotated file in build dir: $build_dir\n\n");
      exit 1;
    }
    if ($self->max_indels_per_file){
      my $var_count = $self->unique_variant_count('-variant_file'=>$tier1_indel_file);
      if ($var_count > $self->max_indels_per_file){
        $self->status_message("Too many variants ... skipping file: $tier1_indel_file");        
      }else{
        push(@file_list, $tier1_indel_file) unless ($var_count > $var_count);
      }
    }else{
      push(@file_list, $tier1_indel_file);
    }
  }

  #Make sure at least one file made it through the filters if not, exit with status 0
  unless (scalar(@file_list) > 0){
    $self->warning_message("All files were skipped because they exceeded the max specified number of SNVs or INDELs - no plots will be generated");
    exit 0;
  }


  #If the user specified to merge duplicate variants, do so now
  my $file_list_string = join(" ", @file_list);
  my $merge_cmd;
  if ($self->collapse_variants){
    $merge_cmd = "cat $file_list_string | sort | uniq > $complete_variants_file";
  }else{
    $merge_cmd = "cat $file_list_string > $complete_variants_file";
  }
  $self->status_message("\n\n$merge_cmd");
  Genome::Sys->shellcmd(cmd => $merge_cmd);
  
  #Create a filtered version of the variants file that removes variants that match the effects type filter
  open (VAR, "$complete_variants_file") || die "\n\nCould not open complete variants file: $complete_variants_file\n\n";
  open (VARF, ">$filtered_variants_file") || die "\n\nCould not open filtered variants file for writing: $filtered_variants_file\n\n";
  my $filter_string = $self->effect_type_filter;
  my $c = 0;
  while(<VAR>){
    chomp($_);
    my @line = split("\t", $_);
    my $chr = $line[0];
    my $start = $line[1];
    my $end = $line[2];
    my $coord = "$chr:$start-$end";
    my $gid = $line[6];
    my $tid = $line[7];
    my $cdna_position = $tid . "_" . $line[14];
    my $aa_effect = $line[15];
    my $aa_position = $tid . "_" . $aa_effect;
    my $aa_id = $aa_position;
    my $chr_position = $tid . "_" . $coord;
    if ($aa_effect =~ /(\d+)/){
      $aa_position = $tid . "_" . $1;
    }
    
    next if ($line[13] =~ /$filter_string/);
    $c++;

    if ($self->max_transcripts){
      next if ($c > $self->max_transcripts);
    }

    print VARF "$_\n";

    if ($var{$tid}){
      my $transcript_variants = $var{$tid}{variants};
      $transcript_variants->{$c}->{var} = $_;
      $transcript_variants->{$c}->{aa_id} = $aa_id;
      $transcript_variants->{$c}->{cdna_position} = $cdna_position;
      $transcript_variants->{$c}->{chr_position} = $chr_position;
      $transcript_variants->{$c}->{aa_position} = $aa_position;
    }else{
      my %tmp;
      $tmp{$c}{var} = $_;
      $tmp{$c}{aa_id} = $aa_id;
      $tmp{$c}{cdna_position} = $cdna_position;
      $tmp{$c}{chr_position} = $chr_position;
      $tmp{$c}{aa_position} = $aa_position;
      $var{$tid}{variants} = \%tmp;
      $var{$tid}{gid} = $gid;
    }
  }
  close(VAR);
  close (VARF);

  return (\%var);
}



#####################################################################################################################
#Import variants from the the cosmic build                                                                          #
#####################################################################################################################
sub import_cosmic_variants{
  my $self = shift;
  my %args = @_;
  my $somatic_variants = $args{'-somatic_variants'};
  my $cosmic_annotation_file = $args{'-cosmic_annotation_file'};
  my $filtered_variants_file = $args{'-filtered_variants_file'};

  my @builds = $self->builds;
  my $outdir = $self->outdir;
  unless ($outdir =~ /\/$/){
    $outdir .= "/";
  }

  my %trans_var;
  my %cosmic_aa_positions;
  my %cosmic_cdna_positions;
  my %cosmic_chr_positions;
  my %cosmic_aa_ids;

  my $transcript_summary_file = $outdir . "cosmic_transcript_somatic_summary.tsv";
  my $variants_intersect_cosmic_file = $outdir . "variants.hq.tier1.v1.annotated.cosmic.tsv";

  #parse through the cosmic variants and grab those lines that correspond to transcripts in the list with somatic variants
  my $header;
  open (COSMIC, "$cosmic_annotation_file") || die "\n\nCould not open cosmic file: $cosmic_annotation_file\n\n";
  while(<COSMIC>){
    chomp($_);
    my @line = split("\t", $_);
    if ($_ =~ /^chromosome/){
      $header = $_;
      next;
    }
    my $chr = $line[0];
    my $start = $line[1];
    my $end = $line[2];
    my $coord = "$chr:$start-$end";
    my $gid = $line[6];
    my $tid = $line[7];
    my $cdna_position = $tid . "_" . $line[14];
    my $aa_effect = $line[15];
    my $aa_position = $tid . "_" . $aa_effect;
    my $aa_id = $aa_position;
    my $chr_position = $tid . "_" . $coord;
    if ($aa_effect =~ /(\d+)/){
      $aa_position = $tid . "_" . $1;
    }

    $cosmic_aa_ids{$aa_id}++;
    $cosmic_aa_positions{$aa_position}++;
    $cosmic_cdna_positions{$cdna_position}++;
    $cosmic_chr_positions{$chr_position}++;

    #Skip unless the cosmic variants on transcripts that were not in the list of somatic mutated transcripts
    next unless $somatic_variants->{$tid};

    #Hash all annotation lines and keep counts for each *transcript*
    if ($trans_var{$tid}){
      my $anno = $trans_var{$tid}{anno};
      $anno->{$_}++;
      $trans_var{$tid}{cosmic_mutation_sample_count}++;
    }else{
      my %tmp;
      $tmp{$_}=1;
      $trans_var{$tid}{anno} = \%tmp;
      $trans_var{$tid}{gene} = $gid;
      $trans_var{$tid}{cosmic_mutation_sample_count}++; #Total observations of any mutation for this transcript
      $trans_var{$tid}{cosmic_mutation_site_count} = 0; #Total distinct mutated sites for this transcript in cosmic
    }
  }
  close (COSMIC);

  #How many variants are there per transcript, what is the max observation count for a single mutation?
  foreach my $tid (sort {$trans_var{$a}->{gene} cmp $trans_var{$b}->{gene}} keys %trans_var){
    my $anno = $trans_var{$tid}{anno};
    my $gid = $trans_var{$tid}{gene};
    my $mutation_count = keys %{$anno};
    my $max_mut_obs = 0;
    $trans_var{$tid}{cosmic_mutation_site_count} = $mutation_count;
    foreach my $mut (keys %{$anno}){
      my $mut_obs = $anno->{$mut};
      $max_mut_obs = $mut_obs if ($mut_obs > $max_mut_obs);
    }
    $trans_var{$tid}{max_mutations_observed} = $max_mut_obs;
  }

  #Print out a transcript level summary of Cosmic variants for transcripts with a somatic variant

  open (OUT, ">$transcript_summary_file") || die "\n\nCould not open transcript summary file for writing: $transcript_summary_file\n\n";
  print OUT "transcript_id\tgene_id\ttotal_sample_mutation_count\tdistinct_mutation_sites\tmax_mutations_observed_at_a_site\n";
  foreach my $tid (sort {$somatic_variants->{$a}->{gid} cmp $somatic_variants->{$b}->{gid}} keys %{$somatic_variants}){
    if (defined($trans_var{$tid})){
      print OUT "$tid\t$trans_var{$tid}{gene}\t$trans_var{$tid}{cosmic_mutation_sample_count}\t$trans_var{$tid}{cosmic_mutation_site_count}\t$trans_var{$tid}{max_mutations_observed}\n";
    }else{
      print OUT "$tid\t$somatic_variants->{$tid}->{gid}\t0\t0\t0\n";
    }
  }
  close (OUT);

  #Print out the original list of somatic variant annotation records with occurences of that exact amino acid change in Cosmic reported
  #By keying on $tid+$aa_change

  open (OUT, ">$variants_intersect_cosmic_file") || die "\n\nCould not open variants intersect cosmic file for writing: $variants_intersect_cosmic_file";
  print OUT "$header\tcosmic_aa_id_match_count\tcosmic_aa_position_match_count\tcosmic_cdna_position_match_count\tcosmic_chr_position_match_count\n";
  foreach my $tid (sort keys %{$somatic_variants}){
    my $transcript_variants = $somatic_variants->{$tid}->{variants};
    foreach my $c (sort {$a <=> $b} keys %{$transcript_variants}){
      my $aa_id = $transcript_variants->{$c}->{aa_id};
      my $aa_position = $transcript_variants->{$c}->{aa_position};
      my $cdna_position = $transcript_variants->{$c}->{cdna_position};
      my $chr_position = $transcript_variants->{$c}->{chr_position};
      my $line = $transcript_variants->{$c}->{var};
      my $cosmic_aa_id_match = 0;
      my $cosmic_aa_position_match = 0;
      my $cosmic_cdna_position_match = 0;
      my $cosmic_chr_position_match = 0;
      if ($cosmic_aa_ids{$aa_id}){
        $cosmic_aa_id_match = $cosmic_aa_ids{$aa_id};
      }
      if ($cosmic_aa_positions{$aa_position}){
        $cosmic_aa_position_match = $cosmic_aa_positions{$aa_position};
      }
      if ($cosmic_cdna_positions{$cdna_position}){
        $cosmic_cdna_position_match = $cosmic_cdna_positions{$cdna_position};
      }
      if ($cosmic_chr_positions{$chr_position}){
        $cosmic_chr_position_match = $cosmic_chr_positions{$chr_position};
      }
      print OUT "$line\t$cosmic_aa_id_match\t$cosmic_aa_position_match\t$cosmic_cdna_position_match\t$cosmic_chr_position_match\n"
    }
  }

  #If there is a mutation in a transcript with more than N variants, rescale the counts to max out at $max_var (or perhap convert to log2 scale?)
  my $t_count = keys %trans_var;
  $self->status_message("Creating filtered Cosmic variants file for $t_count transcripts that had somatic variants");
  open (OUT, ">$filtered_variants_file") || die "\n\nCould not open cosmic variants file for output: $filtered_variants_file\n\n";
  foreach my $tid (sort {$trans_var{$a}->{gene} cmp $trans_var{$b}->{gene}} keys %trans_var){
    my $anno = $trans_var{$tid}{anno};
    my $gid = $trans_var{$tid}{gene};
    my $mutation_count = keys %{$anno};
    my $max_mutation_count = $trans_var{$tid}{max_mutations_observed};
    my $sites = 0;
    foreach my $record (sort {$anno->{$b} <=> $anno->{$a}} keys %{$anno}){
      $sites++;
      last if ($sites >= $self->max_cosmic_sites);

      my $record_count = $anno->{$record};

      #If there is a max mutation count higher than that allowed, all record counts will be rescaled.  Otherwise they will be unchanged
      if ($max_mutation_count > $self->max_count_per_cosmic_site){
        my $x = ($record_count/$max_mutation_count)*$self->max_count_per_cosmic_site;
        my $xf = sprintf("%.0f", $x);
        $xf = 1 if ($xf < 1);
        $record_count = $xf;
      }
      for (my $i = 1; $i <= $record_count; $i++){
        print OUT "$record\n";
      }
    }
  }
  close (OUT);

  return (\%trans_var);
}





#####################################################################################################################
#Draw mutation diagrams - one for somatic variants specified and one for all cosmic variants                        #
#####################################################################################################################
sub draw_mutation_diagrams{
  my $self = shift;
  my %args = @_;
  my $somatic_variants_file = $args{'-somatic_variants_file'};
  my $cosmic_variants_file = $args{'-cosmic_variants_file'};

  my $ab_name = $args{'-annotation_build_name'};

  my $outdir = $self->outdir;
  unless ($outdir =~ /\/$/){
    $outdir .= "/";
  }

  $outdir .= "images/";
  mkdir($outdir);
  
  my $mut_diag_cmd1 = Genome::Model::Tools::Graph::MutationDiagram->create(annotation=>$somatic_variants_file, reference_transcripts=>$ab_name, output_directory=>$outdir, file_suffix=>'_SOMATIC');
  $mut_diag_cmd1->execute();

  my $mut_diag_cmd2 = Genome::Model::Tools::Graph::MutationDiagram->create(annotation=>$cosmic_variants_file, reference_transcripts=>$ab_name, output_directory=>$outdir, file_suffix=>'_COSMIC');
  $mut_diag_cmd2->execute();
  
  return;
}


#####################################################################################################################
#Count unique variants in an annotation file                                                                        #
#####################################################################################################################
sub unique_variant_count{
  my $self = shift;
  my %args = @_;
  my $variant_file = $args{'-variant_file'};  
  my %vars;
  open (VAR, "$variant_file") || die "\n\nCould not open variant file: $variant_file\n\n";
  while(<VAR>){
    chomp($_);
    my @line = split("\t", $_);
    my $coord = "$line[0]:$line[1]-$line[2]";
    $vars{$coord}=1;
  }
  close(VAR);
  my $variant_count = keys %vars;
  return($variant_count);
}


1;


