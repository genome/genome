package Genome::Model::Command::Services::SummarizeOutstandingBuilds;

use strict;
use warnings;
use Genome;
use Date::Calc qw(Delta_DHMS);

class Genome::Model::Command::Services::SummarizeOutstandingBuilds{
    is => 'Command',
    has => [
        take_action => {
            is => 'Boolean',
            is_input => 1,
            default => 0,
            doc => 'Currently disabled.  Takes appropriate action for each group of models or builds if true.  Otherwise only prints out a report',
        },
    ],
};

sub help_brief {
    #TODO: write me
    return "";
}

sub help_synopsis{
    #TODO: write me
    return "";
}

sub help_detail {
    #TODO: write me
    return "";
}

sub execute {
    my $self = shift;

    my $lock_resource = $ENV{GENOME_LOCK_DIR} . '/genome_model_command_services_summarize-outstanding-builds/loader';
    my $lock = Genome::Sys->lock_resource(resource_lock=>$lock_resource, max_try=>1);
    unless ($lock){
        $self->error_message("could not lock, another instance must be running.");
        return;
    }

    $self->initialize_categories;

    my @models = $self->get_models;
    $self->status_message("Done getting models");
    for my $model (@models){
        $self->categorize_models_for_processing($model);
    }
    $self->status_message("Done separating models");

    $self->process_categories;

    Genome::Sys->unlock_resource(resource_lock =>$lock);
 
    return 1;
}


sub get_models {
    my $self = shift;
    return ($self->get_models_with_unsucceeded_builds, $self->get_models_with_no_builds);
}

sub get_models_with_unsucceeded_builds {
    my $self = shift;
    #The events get here is a to try and avoid an n+1 problem with these data structures.  As of early March 2011, it doesn't work, and makes this slow as hell.
    my @events = Genome::Model::Event->get(event_type => 'genome model build', 'event_status ne' => 'Succeeded');
    my @builds = Genome::Model::Build->get('status ne' => 'Succeeded', -hint => ['model', 'the_master_event']);

    my %models;

    for my $build (@builds){
        $models{$build->model->id} = $build->model; 
    }

    return values %models;
}

sub get_models_with_no_builds {
    my $self = shift;
    my %models_without_builds;
    my @models = Genome::Model->get(-hint => ['builds']);

    for my $model (@models){
        my @builds = $model->builds;
        unless(@builds){
            $models_without_builds{$model->id} = $model;
        }
    }

    return values %models_without_builds;
}

sub process_categories {
    my $self = shift;
    my @categories = $self->_category_names;
    for my $category (@categories){
        my $method_name = "_process" . $category;
        $self->$method_name;
    }
}

sub categorize_models_for_processing {
    my ($self, $model) = @_;
    my $latest_build = $model->latest_build;
    if (not $latest_build){
        $self->categorize_model_with_no_builds($model);            
    }elsif($latest_build->status eq 'Scheduled'){
        $self->categorize_model_with_scheduled_build($model, $latest_build);
    }elsif($latest_build->status eq 'Running'){
        $self->categorize_model_with_running_build($model, $latest_build);
    }elsif($latest_build->status eq 'Succeeded'){
        $self->categorize_model_with_succeeded_build($model, $latest_build);
    }elsif($latest_build->status eq 'Abandoned'){
        $self->categorize_model_with_abandoned_build($model, $latest_build);
    }elsif($latest_build->status eq 'Failed' or $latest_build->status eq 'Crashed'){
        $self->categorize_model_with_failed_build($model, $latest_build);
    }elsif($latest_build->status eq 'Preserved'){
        $self->categorize_model_with_preserved_build($model, $latest_build);
    }else{
        $self->warning_message("I dunno what to do with model " . $model->id);
    }
}

sub categorize_model_with_no_builds{
    my $self = shift;
    my $model = shift;
    my $model_age = $self->_age_in_days($model->creation_date);
    if($model_age > 7){
        push(@{$self->_category_hash->{_none_kill_and_email_user}}, $model);
    }elsif($model_age > 3){
        push(@{$self->_category_hash->{_none_email_user}}, $model);
    }
}

sub categorize_model_with_scheduled_build{
    my ($self, $model, $latest_build) = @_;
    # my $latest_build_id = $latest_build->id;
    # my $lsf_job_id = $latest_build->the_master_event->lsf_job_id;
    
    # my $bjobs_output = `bjobs $lsf_job_id`;
    # if(not $bjobs_output or not ($bjobs_output =~ m/$latest_build_id/)){
        # push(@{$self->_category_hash->{_scheduled_incorrect_state}}, $latest_build);
    # }
#TODO: replace this with something that works.  lsf_job_id isn't useful if the build is scheduled
}

sub categorize_model_with_running_build{
    my ($self, $model, $latest_build) = @_;

    my $latest_build_id = $latest_build->id;
    my $job_id;
    my $workflow = $latest_build->newest_workflow_instance;
    $job_id = $workflow->current->dispatch_identifier unless ($job_id || !$workflow);
    $job_id = $latest_build->the_master_event->lsf_job_id unless ($job_id);
    my $bjobs_output = `bjobs $job_id`;
    if(not $bjobs_output or not ($bjobs_output =~ m/$latest_build_id/)){
        push(@{$self->_category_hash->{_running_incorrect_state}}, $latest_build);
        return;
    }

    my ($model_age) = $self->_age_in_days($latest_build->date_scheduled);
    if($model_age >=10){
        push(@{$self->_category_hash->{_running_kill_and_email_user}}, $latest_build);
    }elsif($model_age >= 3){
        push(@{$self->_category_hash->{_running_email_user}}, $latest_build);
    }
}

sub categorize_model_with_succeeded_build{
    my ($self, $model, $latest_build) = @_;
    my @builds = $model->builds;
    for my $build (@builds){
        next if $build == $latest_build;
        if ($build->status eq 'Failed'){
            push(@{$self->_category_hash->{_succeeded_abandon_builds}}, $build);
        }elsif ($build->status eq 'Succeeded'){
            push(@{$self->_category_hash->{_succeeded_eviscerate_builds}}, $build);
        }elsif ($build->status eq 'Running' and $self->_age_in_days($latest_build->date_completed) >= 1){
            push(@{$self->_category_hash->{_succeeded_kill_running_builds}}, $build) if $self->_age_in_days($build->date_scheduled);
        }
    }
}

sub categorize_model_with_abandoned_build{
    my ($self, $model, $latest_build) = @_;
    push(@{$self->_category_hash->{_abandoned_start_new_builds}}, $model);
}

sub categorize_model_with_failed_build{
    my ($self, $model, $latest_build) = @_;
    my @builds = $model->builds;

    if(scalar(@builds) > 2 and ($builds[-2]->status eq 'Failed' or $builds[-2]->status eq 'Crashed')){
        #flag model for manual attention if the last two builds are failures
        push(@{$self->_category_hash->{_failed_flag_for_manual_attention}}, $model);
    }
    else{
        push(@{$self->_category_hash->{_failed_start_new_builds}}, $model);
    }
}

sub categorize_model_with_preserved_build{
    my ($self, $model, $latest_build) = @_;
    #TODO: what do we do?  
    push(@{$self->_category_hash->{_preserved}}, $model);
}

sub initialize_categories {
    my $self = shift;
    my %category_hash;
    for my $category ($self->_category_names){
        $category_hash{$category} = [];
    }
    
    $self->_category_hash(\%category_hash);
    return 1;
}

sub _process_none_email_user{
    my $self = shift;
    if($self->take_action){
        #TODO: iterate the category and take appropriate action 
    }

    print "E-MAIL FOR MODEL WITH NO BUILDS\n";
    for my $model (@{$self->_category_hash->{_none_email_user}}){
        print join("\t", $model->id, $model->creation_date), "\n"; 
    }
}

sub _process_none_kill_and_email_user{
    my $self = shift;
    if($self->take_action){
        #TODO: iterate the category and take appropriate action 
    }

    print "KILL AND E-MAIL FOR MODEL WITH NO BUILDS\n";
    for my $model (@{$self->_category_hash->{_none_kill_and_email_user}}){
        print join("\t", $model->id, $model->creation_date), "\n"; 
    }
}

sub _process_scheduled_incorrect_state {
    my $self = shift;
    if($self->take_action){
        #TODO: iterate the category and take appropriate action 
    }

    print "SCHEDULED BUILDS IN INCORRECT STATE\n";
    for my $build (@{$self->_category_hash->{_scheduled_incorrect_state}}){
        print join("\t", $build->model->id, $build->id, $build->date_scheduled, $build->status), "\n"; 
    }
}

sub _process_running_incorrect_state{
    my $self = shift;
    if($self->take_action){
        #TODO: iterate the category and take appropriate action 
    }

    print "RUNNING BUILDS IN INCORRECT STATE\n";
    for my $build (@{$self->_category_hash->{_running_incorrect_state}}){
        print join("\t", $build->model->id, $build->id, $build->date_scheduled, $build->status), "\n"; 
    }
}

sub _process_running_email_user{
    my $self = shift;
    if($self->take_action){
        #TODO: iterate the category and take appropriate action 
    }

    print "EMAIL ABOUT RUNNING BUILDS 3-9 DAYS OLD\n";
    for my $build (@{$self->_category_hash->{_running_email_user}}){
        print join("\t", $build->model->id, $build->id, $build->date_scheduled, $build->status), "\n"; 
    }
}

sub _process_running_kill_and_email_user{
    my $self = shift;
    if($self->take_action){
        #TODO: iterate the category and take appropriate action 
    }

    print "KILL AND EMAIL ABOUT RUNNING BUILDS 10+ DAYS OLD\n";
    for my $build (@{$self->_category_hash->{_running_kill_and_email_user}}){
        print join("\t", $build->model->id, $build->id, $build->date_scheduled, $build->status), "\n"; 
    }
}

sub _process_succeeded_abandon_builds{
    my $self = shift;
    if($self->take_action){
        #TODO: iterate the category and take appropriate action 
    }

    print "SUCCEEDED BUILD WITH OLD FAILED BUILD TO ABANDON\n";
    for my $build (@{$self->_category_hash->{_succeeded_abandon_builds}}){
        print join("\t", $build->model->id, $build->id, $build->date_scheduled, $build->status), "\n"; 
    }
}

sub _process_succeeded_eviscerate_builds{
    my $self = shift;
    if($self->take_action){
        #TODO: iterate the category and take appropriate action 
    }

    print "SUCCEEDED BUILD WITH OLD SUCCEEDED BUILD TO EVISCERATE\n";
    for my $build (@{$self->_category_hash->{_succeeded_eviscerate_builds}}){
        print join("\t", $build->model->id, $build->id, $build->date_scheduled, $build->status), "\n"; 
    }
}

sub _process_succeeded_kill_running_builds{
    my $self = shift;
    if($self->take_action){
        #TODO: iterate the category and take appropriate action 
    }

    print "SUCCEEDED BUILD WITH OLD RUNNING BUILD TO KILL\n";
    for my $build (@{$self->_category_hash->{_succeeded_kill_running_builds}}){
        print join("\t", $build->model->id, $build->id, $build->date_scheduled, $build->status), "\n"; 
    }
}

sub _process_abandoned_start_new_builds{
    my $self = shift;
    if($self->take_action){
        #TODO: iterate the category and take appropriate action 
    }

    print "START BUILD ON MODEL WITH ABANDONED BUILD\n";
    for my $model (@{$self->_category_hash->{_abandoned_start_new_builds}}){
        print join("\t", $model->id, $model->creation_date), "\n"; 
    }
}

sub _process_failed_start_new_builds{
    my $self = shift;
    if($self->take_action){
        #TODO: iterate the category and take appropriate action 
    }

    print "START BUILD ON MODEL WITH SINGLE FAILED BUILD\n";
    for my $model (@{$self->_category_hash->{_failed_start_new_builds}}){
        print join("\t", $model->id, $model->creation_date), "\n"; 
    }
}

sub _process_failed_flag_for_manual_attention{
    my $self = shift;
    if($self->take_action){
        #TODO: iterate the category and take appropriate action 
    }

    print "FLAG MODEL WITH TWO FAILED BUILDS IN A ROW FOR MANUAL ATTENTION\n";
    for my $model (@{$self->_category_hash->{_failed_flag_for_manual_attention}}){
        print join("\t", $model->id, $model->creation_date), "\n"; 
    }
}

sub _process_preserved {
    my $self = shift;
    if($self->take_action){
        #TODO: iterate the category and take appropriate action 
    }

    print "IGNORE MODEL WITH PRESERVED BUILD\n";
    for my $model (@{$self->_category_hash->{_preserved}}){
        print join("\t", $model->id, $model->creation_date), "\n"; 
    }
}

sub _age_in_days {
    my ($self, $date) = @_;
    my ($age_in_days) = Delta_DHMS(split("-|:| ", $date), split("-|:| ",$self->__context__->now));
    return $age_in_days;
}

sub _category_names{
    my @categories = qw(_none_email_user _none_kill_and_email_user _scheduled_incorrect_state _running_incorrect_state _running_email_user _running_kill_and_email_user _succeeded_abandon_builds _succeeded_eviscerate_builds _succeeded_kill_running_builds _abandoned_start_new_builds _failed_start_new_builds _failed_flag_for_manual_attention _preserved);
    return @categories;
}

sub _category_hash{
    my ($self, $hash_ref) = @_;
    $self->{_category_hash} = $hash_ref if $hash_ref;
    return $self->{_category_hash};
}
