package Genome::Model::ClinSeq::Command::Converge::Base;
use strict;
use warnings;
use Genome;
use Genome::Model::ClinSeq; #Needed temporarily to deal with the hack to add build->common_name to the clinseq model/build object
use Data::Dumper;
use List::MoreUtils qw/ uniq /;

class Genome::Model::ClinSeq::Command::Converge::Base {
    is => 'Command::V2',
    is_abstract => 1,
    has_input => [
        builds => { 
            is => 'Genome::Model::Build::ClinSeq',
            is_many => 1,
            require_user_verify => 0,
            doc => 'clinseq builds to converge',
        },
        bam_readcount_version => {
            is => 'String',
            doc => 'version of bam-readcount to use',
        },
        outdir => {
            is => 'FilesystemPath',
            doc => 'Directory where output files will be written',
        },
        mq => {
              is => 'Integer',
              is_optional => 1,
              doc => 'minimum mapping quality of reads to be considered',
              default => '30',
        },
        bq => {
              is => 'Integer',
              is_optional => 1,
              doc => 'minimum base quality of bases in reads to be considered',
              default => '20',
        },
    ],
    doc => 'converge various data types across clinseq inputs'
};

sub __errors__ {
  my $self = shift;
  my @errors = $self->SUPER::__errors__(@_);

  if ($self->outdir){
    unless (-e $self->outdir && -d $self->outdir) {
      push @errors, UR::Object::Tag->create(
                                            type => 'error',
                                            properties => ['outdir'],
                                            desc => "Outdir: " . $self->outdir . " not found or not a directory",
                                          );
    }
  }
  return @errors;
}


#Resolve human readable subject labels (one label per clinseq build)
#Test clin-seq model-groups 51042, 50714, 62686, 52824
sub resolve_clinseq_subject_labels{
  my $self = shift;
  my @builds = $self->builds;
  
  $self->debug_message("Attempting to resolve distinct clinseq subject names in a human readable format");

  my $build_count = scalar(@builds);

  my %patient_common_names;
  my %patient_and_subject_type_names;
  my %clinseq_subject_names;
  my %patient_names; 
  my %wgs_subject_names; 
  my %exome_subject_names;
  my %dna_subject_names;
  my %tumor_rna_subject_names;
  my %normal_rna_subject_names;
  my %model_id_names;

  my @name_lists;
  foreach my $build (@builds){

    #Patient common names. e.g. PNC6
    my $patient_common_name = $self->get_final_common_name('-clinseq_build'=>$build);
    $patient_common_names{$patient_common_name} = $build->id if $patient_common_name;

    #Patient common name combined with dna sample type.  e.g. PNC6_tumor
    my $dna_subject_type = $self->get_dna_sample_type('-clinseq_build'=>$build);
    my $patient_and_subject_type_name;
    if ($patient_common_name && $dna_subject_type){
      $patient_and_subject_type_name = $patient_common_name . "_" . $dna_subject_type;
    }
    $patient_and_subject_type_names{$patient_and_subject_type_name} = $build->id if $patient_and_subject_type_name;

    #Clinseq subject name. e.g. H_KU-321
    my $clinseq_subject_name = $build->subject->name if $build->subject;
    $clinseq_subject_names{$clinseq_subject_name} = $build->id if $build->subject;

    #WGS subject name. e.g. H_KA-306905-1121472
    my $wgs_subject_name = $build->wgs_build->subject->name if $build->wgs_build;
    $wgs_subject_names{$wgs_subject_name} = $build->id if $build->wgs_build;

    #Exome subject name. e.g. H_KA-306905-1121472
    my $exome_subject_name = $build->exome_build->subject->name if $build->exome_build;
    $exome_subject_names{$exome_subject_name} = $build->id if $build->exome_build;

    #DNA subject name. e.g. H_KA-306905-1121472
    my $dna_subject_name = $self->get_dna_subject_name('-clinseq_build'=>$build, '-patient_common_name'=>$patient_common_name);
    $dna_subject_names{$dna_subject_name} = $build->id if $dna_subject_name;

    #Tumor RNA subject name
    my $tumor_rna_subject_name = $build->tumor_rnaseq_build->subject->name if $build->tumor_rnaseq_build;
    $tumor_rna_subject_names{$tumor_rna_subject_name} = $build->id if $build->tumor_rnaseq_build;

    #Normal RNA subject name
    my $normal_rna_subject_name = $build->normal_rnaseq_build->subject->name if $build->normal_rnaseq_build;
    $normal_rna_subject_names{$normal_rna_subject_name} = $build->id if $build->normal_rnaseq_build;

    #Model ID names
    my $model_id = $build->model->id;
    $model_id_names{$model_id} = $build->id;

  }
  push (@name_lists, \%patient_common_names);
  push (@name_lists, \%patient_and_subject_type_names);
  push (@name_lists, \%clinseq_subject_names);
  push (@name_lists, \%wgs_subject_names);
  push (@name_lists, \%exome_subject_names);
  push (@name_lists, \%dna_subject_names);
  push (@name_lists, \%tumor_rna_subject_names);
  push (@name_lists, \%normal_rna_subject_names);
  push (@name_lists, \%model_id_names);

  my $labels1;

  foreach my $name_list (@name_lists){
    if (scalar(keys %{$name_list}) == $build_count){
      $labels1 = $name_list;
      last;
    }
  }

  #If a unique set of names still could not be found use the build ids as a last resort
  my %tmp;
  unless (scalar(keys %{$labels1}) == $build_count){
    foreach my $build (@builds){
      my $build_id = $build->id;
      $tmp{$build_id} = $build_id;
    }
    $labels1 = \%tmp;
  }

  #reverse the key value so that labels are keyed on build id before returning
  my %labels2;
  my $o = 0;
  foreach my $label (sort keys %{$labels1}){
    $o++;
    my $build_id = $labels1->{$label};
    $labels2{$build_id}{name} = $label;
    $labels2{$build_id}{name_abr} = $label;
    $labels2{$build_id}{order} = $o;
  }

  return(\%labels2);
}


sub get_final_common_name{
  my $self = shift;
  my %args = @_;
  my $clinseq_build = $args{'-clinseq_build'};

  my $final_name;
  
  my $wgs_build = $clinseq_build->wgs_build;
  my $exome_build = $clinseq_build->exome_build;
  my $normal_rnaseq_build = $clinseq_build->normal_rnaseq_build;
  my $tumor_rnaseq_build = $clinseq_build->tumor_rnaseq_build;

  my @builds = ($clinseq_build, $wgs_build, $exome_build, $normal_rnaseq_build, $tumor_rnaseq_build);

  my %names;
  foreach my $build (@builds){
    next unless $build;
    if ($build->subject->class eq 'Genome::Individual'){
      my $common_name = $build->subject->common_name;
      $names{$common_name}=1 if $common_name;
      $final_name = $common_name if $common_name;
    }elsif ($build->subject->class eq 'Genome::Sample'){
      my $common_name = $build->subject->individual->common_name;
      $names{$common_name}=1 if $common_name;
      $final_name = $common_name if $common_name;
    }
  }
  if (scalar (keys %names) > 1){
    $self->warning_message("Found multiple patient common names for clin-seq build: " . $clinseq_build->id);
  }

  return ($final_name);
}


sub get_final_name{
  my $self = shift;
  my %args = @_;
  my $clinseq_build = $args{'-clinseq_build'};

  my $final_name;
  
  my $wgs_build = $clinseq_build->wgs_build;
  my $exome_build = $clinseq_build->exome_build;
  my $normal_rnaseq_build = $clinseq_build->normal_rnaseq_build;
  my $tumor_rnaseq_build = $clinseq_build->tumor_rnaseq_build;

  my @builds = ($clinseq_build, $wgs_build, $exome_build, $normal_rnaseq_build, $tumor_rnaseq_build);

  my %names;
  foreach my $build (@builds){
    next unless $build;
    if ($build->subject->class eq 'Genome::Individual'){
      my $name = $build->subject->name;
      $names{$name}=1 if $name;
      $final_name = $name if $name;
    }elsif ($build->subject->class eq 'Genome::Sample'){
      my $name = $build->subject->individual->name;
      $names{$name}=1 if $name;
      $final_name = $name if $name;
    }
  }
  if (scalar (keys %names) > 1){
    $self->warning_message("Found multiple patient names for clin-seq build: " . $clinseq_build->id);
  }

  return ($final_name);
}


sub get_dna_sample_type{
  my $self = shift;
  my %args= @_;
  my $clinseq_build = $args{'-clinseq_build'};

  my $sample_type;
  my $wgs_build = $clinseq_build->wgs_build;
  my $exome_build = $clinseq_build->exome_build;
  my @builds = ($wgs_build, $exome_build);

  my $final_name;
  my %names;
  foreach my $build (@builds){
    next unless $build;
    if ($build->subject->class eq 'Genome::Sample'){
      my $sample_common_name = $build->subject->common_name;
      $names{$sample_common_name}=1 if $sample_common_name;
      $final_name = $sample_common_name if $sample_common_name;
    }
  }
  if (scalar (keys %names) > 1){
    $self->warning_message("Found multiple sample common names for dna subjects of somatic builds for clinseq build: " . $clinseq_build->id);
  }

  $final_name =~ s/ /_/g;

  return ($final_name);
}


sub get_dna_subject_name{
  my $self = shift;
  my %args = @_;
  my $clinseq_build = $args{'-clinseq_build'};
  my $patient_common_name = $args{'-patient_common_name'};

  my $subject_name;
  my $wgs_build = $clinseq_build->wgs_build;
  my $exome_build = $clinseq_build->exome_build;
  my @builds = ($wgs_build, $exome_build);

  my $final_name;
  my %names;
  foreach my $build (@builds){
    next unless $build;
    if ($build->subject->class eq 'Genome::Sample'){
      my $sample_name = $build->subject->name;
      my $sample_common_name = $build->subject->common_name;
      if ($sample_name && $patient_common_name && $sample_common_name){
        $sample_name = $patient_common_name . "_" . $sample_name . "_" . $sample_common_name;
      }elsif($sample_name && $sample_common_name){
        $sample_name = $sample_name . "_" . $sample_common_name;
      }elsif ($sample_name && $patient_common_name){
        $sample_name = $patient_common_name . "_" . $sample_name;
      }
      $names{$sample_name}=1 if $sample_name;
      $final_name = $sample_name if $sample_name;
    }
  }
  if (scalar (keys %names) > 1){
    $self->warning_message("Found multiple sample names for dna subjects of somatic builds for clinseq build: " . $clinseq_build->id);
  }

  $final_name =~ s/ /_/g;

  return ($final_name);
}


sub resolve_input_builds{
  my $self = shift;
  my %args = @_;
  
  #Get all underlying builds that could make up a single clinseq model/build 
  my $clinseq_build = $args{'-clinseq_build'};
  my $wgs_build = $clinseq_build->wgs_build;
  my $exome_build = $clinseq_build->exome_build;
  my $tumor_rnaseq_build = $clinseq_build->tumor_rnaseq_build;
  my $normal_rnaseq_build = $clinseq_build->normal_rnaseq_build;
  my $wgs_normal_refalign_build = $wgs_build->normal_build if ($wgs_build);
  my $wgs_tumor_refalign_build = $wgs_build->tumor_build if ($wgs_build);
  my $exome_normal_refalign_build = $exome_build->normal_build if ($exome_build);
  my $exome_tumor_refalign_build = $exome_build->tumor_build if ($exome_build);
  my @builds = ($wgs_build, $exome_build, $tumor_rnaseq_build, $normal_rnaseq_build, $wgs_normal_refalign_build, $wgs_tumor_refalign_build, $exome_normal_refalign_build, $exome_tumor_refalign_build);
  my @defined_builds;
  foreach my $build (@builds){
    next unless $build;
    push (@defined_builds, $build);
  }
  return(\@defined_builds);
}

sub resolve_clinseq_reference_build{
  my $self = shift;
  my @clinseq_builds = $self->builds;

  $self->debug_message("Attempting to resolve a distinct reference sequence build from the input models of all clinseq builds");
  my $reference_build;
  my %reference_builds;

  foreach my $clinseq_build (@clinseq_builds){
    my @builds = @{$self->resolve_input_builds('-clinseq_build'=>$clinseq_build)};
    foreach my $build (@builds){
      my $m = $build->model;
      if ($m->can("reference_sequence_build")){
        $reference_build = $m->reference_sequence_build;
        my $rb_id = $reference_build->id;
        $reference_builds{$rb_id}=1;
      }
    }
  }
  my $rb_count = keys %reference_builds;
  unless ($rb_count == 1){
    print Dumper %reference_builds;
    die $self->error_message("Found $rb_count reference builds for this group of input builds - must be only one");
  }
  $self->debug_message("Found 1: " . $reference_build->__display_name__);

  return ($reference_build);
}


sub resolve_clinseq_annotation_build{
  my $self = shift;
  my @clinseq_builds = $self->builds;

  $self->debug_message("Attempting to resolve a distinct reference annotation build from the input models of all clinseq builds");
  my $annotation_build;
  my %annotation_builds;

  foreach my $clinseq_build (@clinseq_builds){
    my @builds = @{$self->resolve_input_builds('-clinseq_build'=>$clinseq_build)};
    foreach my $build (@builds){
      my $m = $build->model;
      if ($m->can("annotation_build")){
        $annotation_build = $m->annotation_build;
        my $ab_id = $annotation_build->id;
        $annotation_builds{$ab_id}=1;
      }
      if ($m->can("annotation_reference_build")){
        $annotation_build = $m->annotation_reference_build;
        my $ab_id = $annotation_build->id;
        $annotation_builds{$ab_id}=1;
      }
    }
  }
  my $ab_count = keys %annotation_builds;
  unless ($ab_count == 1){
    print Dumper %annotation_builds;
    die $self->error_message("Found $ab_count annotation builds for this group of input builds - must be only one");
  }
  $self->debug_message("Found 1: " . $annotation_build->__display_name__ . " (" . $annotation_build->name . ")");

  return ($annotation_build);
}


sub get_clinseq_files{
  my $self = shift;
  my @builds = $self->builds;
  my %args = @_;
  my $target = $args{'-target'};

  my %files;
  foreach my $build (@builds){
    my $build_id = $build->id;
    my $build_dir = $build->data_directory;
    my $common_name = $build->common_name;
    my $path = $build_dir . "/" . $common_name . "/" . $target;
    die $self->error_message("Could not find expected file: $path") unless (-e $path);
    $files{$build_id}{path} = $path;

    #Get column positions from the header line of this file
  }
  return(\%files);
}

sub get_clinseq_file{
  my $self = shift;
  my %args = @_;
  my $build = $args{'-clinseq_build'};
  my $target = $args{'-target'};

  my %file;
  my $build_dir = $build->data_directory;
  my $common_name = $build->common_name;
  my $path = $build_dir . "/" . $common_name . "/" . $target;
  die $self->error_message("Could not find expected file: $path") unless (-e $path);
  $file{path} = $path;

  #Get column positions from the header line of this file
  my $header = 1;
  my %cols;
  open (IN, $path) || die $self->error_message("Could not open file: $path");
  while(<IN>){
    chomp($_);
    my @line = split("\t", $_);
    if ($header){
      $header = 0;
      my $p = 0;
      foreach my $col (@line){
        $cols{$col}{p} = $p;
        $p++;
      }
      last;
    }
  }
  close(IN);
  $file{columns} = \%cols;

  return(\%file);
}


sub getModelsBuilds{
  my $self = shift;
  my %args = @_;
  my $builds_ref = $args{'-builds'};
  my $models_ref = $args{'-models'};
  my $model_group_id = $args{'-model_group_id'};
  my $partial = $args{'-partial'};
  my $verbose = $args{'-verbose'};

  if (defined($partial)){
    $partial = 1;
  }else{
    $partial = 0;
  }

  #The user must specify only one of these options for simplicity
  my $opt_count = 0;
  if (defined($builds_ref)){$opt_count++;}
  if (defined($models_ref)){$opt_count++;}
  if (defined($model_group_id)){$opt_count++;}
  if ($opt_count == 0){die $self->error_message("You must specify at least one option for &getModels (Coverge.pm)"); }
  if ($opt_count > 1){die $self->error_message("You must specify only one option for &getModels (Coverge.pm)"); }

  #Determine the target number of models/builds
  my $target_count = 0;
  my $mg;
  if (defined($builds_ref)){$target_count = scalar(@$builds_ref);}
  if (defined($models_ref)){$target_count = scalar(@$models_ref);}
  if (defined($model_group_id)){
    $mg = Genome::ModelGroup->get("id"=>$model_group_id);
    my @models = $mg->models;
    $target_count = scalar(@models);
  }
  if ($verbose){$self->debug_message("Searching for $target_count models/builds");}

  #Always returns model AND build objects for convenience
  #If a model-group ID or array of model IDs is provided, get the last successful build for each and return that
  #If an array of build IDs is provided, get the model from these builds and return the build objects specified (regardless of whether they are the latest)
  #Make sure all builds are successful.  If a model does not have at least once successful model, warn the user
  my @models;
  my @builds;

  #Get models and builds starting with a list of model build IDs
  if (defined($builds_ref)){
    foreach my $build_id (@{$builds_ref}){
      my $b = Genome::Model::Build->get("id"=>$build_id);
      my $m = $b->model;

      #Make sure the build is successful
      my $status = $b->status;
      if ($status eq "Succeeded"){
        push(@builds, $b);
        push(@models, $m);
      }else{
        $self->debug_message("Warning - build $build_id has a status of $status");
      }
    }
  }

  if (defined($models_ref)){
    foreach my $model_id (@{$models_ref}){
      my $m = Genome::Model->get("id"=>$model_id);
      my $b = $m->last_complete_build;
      if ($b){
        my $status = $b->status;
        my $build_id = $b->id;
        if ($status eq "Succeeded"){
          push(@builds, $b);
          push(@models, $m);
        }else{
          $self->warning_message("Warning - build $build_id has a status of $status");
        }
      }else{
        $self->warning_message("Warning - model $model_id has no complete builds");
      }
    }
  }

  if (defined($model_group_id)){
    @models = $mg->models;
    foreach my $m (@models){
      my $model_id = $m->id;
      my $b = $m->last_complete_build;
      if ($b){
        my $status = $b->status;
        my $build_id = $b->id;
        if ($status eq "Succeeded"){
          push(@builds, $b);
        }else{
          $self->warning_message("Warning - build $build_id has a status of $status");
        }
      }else{
        $self->warning_message("Warning - model $model_id has no complete builds");
      }
    }
  }

  #Check processing profile used for all models and warn if the models are not consistent
  my %pp_list;
  foreach my $m (@models){
    my $pp_id = $m->processing_profile_id;
    $pp_list{$pp_id}=1;
  }
  my $pp_count = keys %pp_list;
  if ($pp_count > 1){
    my @pp = keys %pp_list;
    my $pp_list = join(",", @pp);
    $self->warning_message("Warning - more than one processing-profile is being used by this group of models. PP list: $pp_list");
  }

  #Summarize builds found
  my $b_count = scalar(@builds);
  my $m_count = scalar(@models);

  #Allow the user to return an incomplete list, but only if they specify this option...
  unless ($partial){
    unless ($b_count == $target_count && $m_count == $target_count){
      die $self->error_message("Did not find the correct number of successful models/builds");
    }
  }
  if ($verbose){$self->debug_message("\n\tFound $b_count builds and $m_count models");}

  #Build a hash that groups each model/build pair together
  my %mb;
  for (my $i = 1; $i <= $b_count; $i++){
    $mb{$i}{model} = $models[$i-1];
    $mb{$i}{build} = $builds[$i-1];
  }

  my %models_builds;
  $models_builds{models} = \@models;
  $models_builds{builds} = \@builds;
  $models_builds{cases} = \%mb;

  return(\%models_builds);
}


sub get_case_name{
  my $self = shift;
  my @builds = $self->builds;

  #First attempt to find a common name
  my %common_names;
  my $final_common_name;
  foreach my $build (@builds){
    my $name = $self->get_final_common_name('-clinseq_build'=>$build);
    $common_names{$name}=1 if $name;
    $final_common_name = $name;
  }
  my $ncn = keys %common_names;
  if ($ncn > 1){
    die $self->error_message("$ncn cases found among these builds, this tool is meant to operate on builds from a single individual");
  }

  #Second attempt to find a name
  my %names;
  my $final_name;
  foreach my $build (@builds){
    my $name = $self->get_final_name('-clinseq_build'=>$build);
    $names{$name}=1 if $name;
    $final_name = $name;
  }
  my $nn = keys %names;
  if ($nn > 1){
    die $self->error_message("$nn cases found among these builds, this tool is meant to operate on builds from a single individual");
  }

  my $resolved_name;
  if ($final_common_name){
    $resolved_name = $final_common_name;
  }elsif($final_name){
    $resolved_name = $final_name
  }else{
    die $self->error_message("could not find an individual common_name or name in these builds");
  }

  return $resolved_name;
}

sub fill_missing_output {
  my $self = shift;
  my $output_file = shift;
  #If there are missing cells relative to the header, fill in with 'NA's
  my $tmp_file = $output_file . ".tmp";
  open (TMP_IN, "$output_file") || die $self->error_message("Could not open file: $output_file");
  open (TMP_OUT, ">$tmp_file" ) || die $self->error_message("Could not open file: $tmp_file");
  my $header = 1;
  my $target_cols;
  while(<TMP_IN>){
    chomp $_;
    my @line = split("\t", $_);
    if ($header){
      $target_cols = scalar(@line);
      $header = 0;
      print TMP_OUT "$_\n";
      next;
    }
    if (scalar(@line) == $target_cols){
      print TMP_OUT "$_\n";
    }elsif(scalar(@line) < $target_cols){
      my $diff = $target_cols - (scalar(@line));
      my @newvals;
      push @newvals, 'NA' for (1..$diff);
      my $newvals_string = join("\t", @newvals);
      print TMP_OUT "$_\t$newvals_string\n";
    }else{
      die $self->error_message("File has more data value in a row than names in the header");
    }
  }
  close(TMP_IN);
  close(TMP_OUT);
  my $mv_cmd = "mv $tmp_file $output_file";
  Genome::Sys->shellcmd(cmd => $mv_cmd);
}

sub add_read_counts{
  my $self = shift;
  my %args = @_;
  my $align_builds = $args{'-align_builds'};
  my $grand_anno_file = $args{'-anno_file'};
  my $indel_size_limit = 25; #max size of indels to report counts for.
  my ($b_quality, $m_quality);
  $m_quality = $self->mq;
  $b_quality = $self->bq;
  my (@prefixes, @bam_files);
  foreach my $name (sort {$align_builds->{$a}->{order} <=> $align_builds->{$b}->{order}} keys %{$align_builds}){
    push(@prefixes, $align_builds->{$name}->{prefix});
    push(@bam_files, $align_builds->{$name}->{bam_path});
  }
  my $header_prefixes = join(",", @prefixes);
  my $bam_list = join(",", @bam_files);

  #Get the reference fasta
  my $reference_build = $self->resolve_clinseq_reference_build;
  my $reference_fasta = $reference_build->full_consensus_path('fa');

  #gmt analysis coverage add-readcounts --bam-files=? --genome-build=? --output-file=? --variant-file=? [--header-prefixes=?] 
  my $output_file = $self->outdir . "variants.all.anno.readcounts";
  if (-e $output_file){
    $self->warning_message("using pre-generated bam read count file: $output_file");
    return($output_file)
  }

  $self->status_message("genome_build: $reference_fasta");
  $self->status_message("indel_size_limit: $indel_size_limit");
  $self->status_message("min_quality_score: $m_quality");
  $self->status_message("min_base_quality: $b_quality");

  my $add_count_cmd = Genome::Model::Tools::Analysis::Coverage::AddReadcounts->create(
    bam_files=>$bam_list,
    genome_build=>$reference_fasta,
    output_file=>$output_file,
    variant_file=>$grand_anno_file,
    header_prefixes=>$header_prefixes,
    indel_size_limit => $indel_size_limit,
    min_quality_score => $m_quality,
    min_base_quality => $b_quality,
    bam_readcount_version => $self->bam_readcount_version,
  );
  my $r = $add_count_cmd->execute();
  die $self->error_message("add-readcounts cmd unsuccessful") unless ($r);
  $self->fill_missing_output($output_file);
  return ($output_file);
}

#Parse known druggable files
sub parseKnownDruggableFiles{
  my $self = shift;
  my %args = @_;
  my $files = $args{'-files'};
  my @event_types = @{$args{'-event_types'}};
  $self->status_message("Parsing files containing known drug-gene interaction data");

  #Store all results organized by gene and separately by drug-gene interaction and separately by patient
  my %result;
  my %genes;
  my %interactions;
  my %patients;
  my %data_type_sum;

  #Note that the druggable files are already filtered down to only the variant affected genes with a drug interaction
  #To get a sense of the total number of events will have to wait until the annotation files are being proccessed

  foreach my $patient (keys %{$files}){
    $self->status_message("$patient");
    foreach my $event_type (@event_types){
      foreach my $data_type (sort keys %{$files->{$patient}->{$event_type}}){
        my $drug_file_path = $files->{$patient}->{$event_type}->{$data_type}->{drug_file_path};
        $self->status_message("$drug_file_path");
        open (IN, "$drug_file_path") || die "Could not open gene-drug interaction file: $drug_file_path";
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
          my $gene_name = $line[$columns{'gene_name'}{position}];
          my $drug_name = $line[$columns{'drug_name'}{position}];
          my $interaction = "$gene_name"."_"."$drug_name";
          $interactions{$interaction}{gene_name} = $gene_name;
          $interactions{$interaction}{drug_name} = $drug_name;
          my $drug_class = "unknown";
          if (defined($columns{'drug_class'})){
            $drug_class = $line[$columns{'drug_class'}{position}];
          }

          #Store drug class list from santa monica db
          if (defined($genes{$gene_name}{drug_class})){
            my $classes = $genes{$gene_name}{drug_class};
            $classes->{$drug_class} = 1;
          }else{
            my %classes;
            $classes{$drug_class} = 1;
            $genes{$gene_name}{drug_class} = \%classes;
          }

          #If the gene has any events it will be associated with all drugs that interact with that gene
          if (defined($genes{$gene_name}{drug_list})){
            my $drugs = $genes{$gene_name}{drug_list};
            $drugs->{$drug_name} = 1;
          }else{
            my %drugs;
            $drugs{$drug_name} = 1;
            $genes{$gene_name}{drug_list} = \%drugs;
          }

          #Store unique combinations of patient, event type, data type and gene (e.g., PNC4  snv  wgs  KRAS)
          $data_type_sum{$patient}{$event_type}{$data_type}{$gene_name} = 1;

          #Add patient lists specific to this event type
          if (defined($genes{$gene_name}{$event_type})){
            my $patients = $genes{$gene_name}{$event_type}{patient_list};
            $patients->{$patient} = 1;
          }else{
            my %patients;
            $patients{$patient} = 1;
            $genes{$gene_name}{$event_type}{patient_list} = \%patients;
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
          if (defined($genes{$gene_name}{grand_list})){
            my $patients = $genes{$gene_name}{grand_list};
            $patients->{$patient} = 1;
          }else{
            my %patients;
            $patients{$patient} = 1;
            $genes{$gene_name}{grand_list} = \%patients;
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
  $result{genes} = \%genes;
  $result{interactions} = \%interactions;
  $result{data_type_sum} = \%data_type_sum;
  return(\%result);
}

sub get_somatic_subject_common_name {
  my $self = shift;
  my $b = shift;
  my $subject_common_name = "null";
  if($b->model->wgs_model) {
    $subject_common_name = $b->model->wgs_model->last_succeeded_build->subject->common_name;
  } elsif($b->model->exome_model) {
    my $subject_common_name = $b->model->exome_model->last_succeeded_build->subject->common_name;
  }
  return $subject_common_name;
}

#Get input files to be parsed
sub getFiles{
  my $self = shift;
  my %args = @_;
  my $builds = $args{'-builds'};
  my %files;

  $self->status_message("Get annotation files and drug-gene interaction files from these builds");
  foreach my $b (@$builds){
    my $build_directory = $b->data_directory;
    my $subject_common_name = $b->subject->common_name;
    my $subject_name = $self->get_somatic_subject_common_name($b);
    $subject_name =~ s/[\s-]/_/g;
    my $build_id = $b->id;

    #If the subject name is not defined, die
    unless ($subject_name){
      die $self->error_message("Could not determine subject name for build: $build_id");
    }

    my $final_name = "Unknown";
    if ($subject_name) {
      $final_name = $subject_name;
    }
    if ($subject_common_name and $subject_name) {
      $final_name = $subject_common_name . "_" . $subject_name;
    }
    $self->status_message("$final_name\t$build_id\t$build_directory");

    #Some event types could have come from exome, wgs, or wgs_exome... depending on the event type allow these options and check in order
    #1.) Look for SNV files
    my $exome_snv_dgidb = $b->exome_snv_dgidb_file;
    my $wgs_snv_dgidb = $b->wgs_snv_dgidb_file;
    my $exome_snv_annot = $b->exome_snv_tier1_annotated_compact_catanno_file;
    my $wgs_snv_annot = $b->wgs_snv_tier1_annotated_compact_catanno_file;
    if (-e $wgs_snv_dgidb and -e $wgs_snv_annot) {
      $files{$final_name}{snv}{wgs}{drug_file_path} =
        $wgs_snv_dgidb;
      $files{$final_name}{snv}{wgs}{annot_file_path} =
        $wgs_snv_annot;
    }
    if (-e $exome_snv_dgidb and -e $exome_snv_annot){
      $files{$final_name}{snv}{exome}{drug_file_path} =
        $exome_snv_dgidb;
      $files{$final_name}{snv}{exome}{annot_file_path} =
        $exome_snv_annot;
    } elsif (not (-e $exome_snv_dgidb and -e $wgs_snv_dgidb)) {
      $self->warning_message("Could not find snv drug-gene file.");
    }

    #2.) Look for InDel files
    my $exome_indel_dgidb = $b->exome_indel_dgidb_file;
    my $wgs_indel_dgidb = $b->wgs_indel_dgidb_file;
    my $exome_indel_annot = $b->exome_indel_tier1_annotated_compact_catanno_file;
    my $wgs_indel_annot = $b->wgs_indel_tier1_annotated_compact_catanno_file;
    if (-e $exome_indel_dgidb and -e $exome_indel_annot){
      $files{$final_name}{indel}{exome}{drug_file_path} =
        $exome_indel_dgidb;
      $files{$final_name}{indel}{exome}{annot_file_path} =
        $exome_indel_annot;
    }
    if (-e $wgs_indel_dgidb and -e $wgs_indel_annot) {
      $files{$final_name}{indel}{wgs}{drug_file_path} =
        $wgs_indel_dgidb;
      $files{$final_name}{indel}{wgs}{annot_file_path} =
        $wgs_indel_annot;
    } elsif (not (-e $exome_indel_dgidb and -e $wgs_indel_dgidb)) {
      $self->warning_message("Could not find indel drug-gene file.");
    }

    #3.) Look for CNV gain files
    my $wgs_cnv_dgidb = $b->wgs_cnv_dgidb_file;
    my $wgs_cnv_annot = $b->wgs_cnv_annot_file;
    if (-e $wgs_cnv_annot and -e $wgs_cnv_dgidb){
      $files{$final_name}{cnv_gain}{wgs}{drug_file_path} =
        $wgs_cnv_dgidb;
      $files{$final_name}{cnv_gain}{wgs}{annot_file_path} =
        $wgs_cnv_annot;
    } else{
      $self->warning_message("Could not find CNV drug-gene file for "
        . "$final_name ($subject_name - $subject_common_name) in:" .
        "\n\t$build_directory");
    }

    #4.) Look for Cufflinks RNAseq outlier expression files
    my $rna_cufflinks_dgidb = $b->rnaseq_tumor_cufflinks_dgidb_file;
    my $rna_cufflinks_annot = $b->rnaseq_tumor_cufflinks_annot_file;
    if (-e $rna_cufflinks_annot and -e $rna_cufflinks_dgidb){
      $files{$final_name}{rna_cufflinks_absolute}{rnaseq}{annot_file_path} =
        $rna_cufflinks_annot;
      $files{$final_name}{rna_cufflinks_absolute}{rnaseq}{drug_file_path} =
        $rna_cufflinks_dgidb;
    }else{
      $self->warning_message("Could not find Cufflinks drug-gene file");
    }

    #5.) Look for Tophat junction RNAseq outlier expression files
    my $rna_tophat_dgidb = $b->rnaseq_tumor_tophat_dgidb_file;
    my $rna_tophat_annot = $b->rnaseq_tumor_tophat_annot_file;
    if (-e $rna_tophat_dgidb and $rna_tophat_annot) {
      $files{$final_name}{rna_tophat_absolute}{rnaseq}{drug_file_path} =
        $rna_tophat_dgidb;
      $files{$final_name}{rna_tophat_absolute}{rnaseq}{annot_file_path} =
        $rna_tophat_annot;
    }else{
      $self->warning_message("Could not find Tophat drug-gene file");
    }
  }

  return(\%files);
}

1;

