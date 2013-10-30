package Genome::Model::ClinSeq::Converge;

#Written by Malachi Griffith

require Exporter;

@ISA = qw( Exporter );
@EXPORT = qw();

@EXPORT_OK = qw(
                &getModelsBuilds
               );

%EXPORT_TAGS = (
                all => [qw(&getModelsBuilds)]
               );

use strict;
use warnings;
use Genome;
use Data::Dumper;
use Genome::Model::ClinSeq::Util qw(:all);


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

