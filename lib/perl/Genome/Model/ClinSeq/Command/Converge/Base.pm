package Genome::Model::ClinSeq::Command::Converge::Base;
use strict;
use warnings;
use Genome;
use Genome::Model::ClinSeq; #Needed temporarily to deal with the hack to add build->common_name to the clinseq model/build object
use Data::Dumper;

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
        outdir => {
                    is => 'FilesystemPath',
                    doc => 'Directory where output files will be written',
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
  
  $self->status_message("Attempting to resolve distinct clinseq subject names in a human readable format");

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
      my $common_name = $build->subject->patient->common_name;
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
      my $name = $build->subject->patient->name;
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


sub resolve_somatic_builds{
  my $self = shift;
  my @clinseq_builds = $self->builds;

  my %somatic_builds;
  foreach my $clinseq_build (@clinseq_builds){
    my $wgs_build = $clinseq_build->wgs_build;
    $somatic_builds{$wgs_build->id}{build} = $wgs_build if $wgs_build;
    $somatic_builds{$wgs_build->id}{type} = 'wgs' if $wgs_build;
    my $exome_build = $clinseq_build->exome_build;
    $somatic_builds{$exome_build->id}{build} = $exome_build if $exome_build;
    $somatic_builds{$exome_build->id}{type} = 'exome' if $exome_build;
  }
  return (\%somatic_builds);
}


sub resolve_clinseq_reference_build{
  my $self = shift;
  my @clinseq_builds = $self->builds;

  $self->status_message("Attempting to resolve a distinct reference sequence build from the input models of all clinseq builds");
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
  $self->status_message("\tFound 1: " . $reference_build->__display_name__);

  return ($reference_build);
}


sub resolve_clinseq_annotation_build{
  my $self = shift;
  my @clinseq_builds = $self->builds;

  $self->status_message("Attempting to resolve a distinct reference annotation build from the input models of all clinseq builds");
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
  $self->status_message("\tFound 1: " . $annotation_build->__display_name__ . " (" . $annotation_build->name . ")");

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

####################################################################################################################
#getModels - get an array of models and builds from a model group id, array of model ids, or array of build ids    #
####################################################################################################################
sub getModelsBuilds{
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
  if ($opt_count == 0){die Genome::Model->error_message("\n\nYou must specify at least one option for &getModels (Coverge.pm)\n\n"); }
  if ($opt_count > 1){die Genome::Model->error_message("\n\nYou must specify only one option for &getModels (Coverge.pm)\n\n"); }

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
  if ($verbose){Genome::Model->status_message("\nSearching for $target_count models/builds");}

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
        Genome::Model->status_message("\n\tWarning - build $build_id has a status of $status");
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
          Genome::Model->warning_message("\n\tWarning - build $build_id has a status of $status");
        }
      }else{
        Genome::Model->warning_message("\n\tWarning - model $model_id has no complete builds");
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
          Genome::Model->warning_message("\n\tWarning - build $build_id has a status of $status");
        }
      }else{
        Genome::Model->warning_message("\n\tWarning - model $model_id has no complete builds");
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
    Genome::Model->warning_message("\n\nWarning - more than one processing-profile is being used by this group of models. PP list: $pp_list");
  }

  #Summarize builds found
  my $b_count = scalar(@builds);
  my $m_count = scalar(@models);

  #Allow the user to return an incomplete list, but only if they specify this option...
  unless ($partial){
    unless ($b_count == $target_count && $m_count == $target_count){
      die Genome::Model->error_message("\n\tDid not find the correct number of successful models/builds");
    }
  }
  if ($verbose){Genome::Model->status_message("\n\tFound $b_count builds and $m_count models");}

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

1;

