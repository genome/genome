package Genome::ModelGroup::Command::SimpleCopy;

use strict;
use warnings;

use Genome;

class Genome::ModelGroup::Command::SimpleCopy {
    is => 'Genome::Command::Base',
    has => [
        from => {
            shell_args_position => 1,
            is => 'Genome::ModelGroup',
            doc => 'the existing group to copy', 
        },
        to => {
            shell_args_position => 2,
            is => 'Text',
            doc => 'the name of the group to create',
        },
        new_pp_id => {
            shell_args_position => 3,
            is => 'Text',
            doc => 'the id of the pp you want to use',
        },
        text_to_append => {
            shell_args_position => 4,
            is => 'Text', 
            doc => "Text to append to all the previous model names. eg if you are copying a model group containing a model named 'PRC1', and supply '_copy' you will get 'PRC1_copy'",
            is_optional=>0,
        }
    ],
    has_optional => {
         },
    doc => 'make a new model group from another, varying properties on the model as specified' 
};

sub help_synopsis {
    return <<EOS
genome model-group copy oldgroup newgroup new_pp_id 

EOS
}


sub execute {
    my $self = shift;

    my $from = $self->from;
    my $to_name = $self->to;
    my $to = Genome::ModelGroup->get(name => $to_name);
    if ($to) {
        die $self->error_message("model group $to_name exists!: " . $to->__display_name__);
    }
    my $pp_id = $self->new_pp_id;
    my $pp = Genome::ProcessingProfile->get($pp_id);
    unless($pp) {
        die $self->error_message("unable to ->get processing profile with id $pp_id");
    }

    Genome::Model->get(id=>[map {$_->model_id} $from->model_bridges]);

    my @from_models = $from->models;

    # Now grab all the instrument data.
    $self->status_message("Preloading instrument data.  This might take a few moments.");
    my @ids = map { $_->instrument_data } @from_models;
  
    # Now preload all the subjects. 
    my %subjects;
    for (@from_models) {
        push @{$subjects{$_->subject_class_name}}, $_->subject_id;
    }

    for my $i (keys %subjects) {
        $self->status_message("Preloading subjects of type $i.");
        my @records = $i->get(id=>$subjects{$i});
    
        # samples get fetched by name, too, so preload them that way
        if ($i eq "Genome::Sample") {
            $self->status_message("Preloading all the samples by name.");
            Genome::Sample->get(name=>[map {$_->name} @records])
        }
    }

    # preload all the model inputs and PP params.
    my @inputs = Genome::Model::Input->get(model_id=>[map {$_->id} @from_models]);
    Genome::ProcessingProfile::Param->get(processing_profile_id=>[map {$_->processing_profile_id} @from_models]);

    # End preloading.  From here the only queries hitting the DB should be to get sequences,
    # and to check the DB for models that already exist with our generated names from this copy.


    $to = Genome::ModelGroup->create(name => $to_name);

    my @new_models;
    my $n = 0;
    my $text_to_append = $self->text_to_append;
    for my $from_model (@from_models) {
        my $from_profile = $from_model->processing_profile;
        my $to_profile = Genome::ProcessingProfile->get($pp_id);
        my $to_model;
            my $new_name;
                $new_name = $from_model->name . "$text_to_append";
             $to_model = $from_model->copy(
                name => $new_name,
                processing_profile => $to_profile
             );            
            
            if($to_model) {
                $self->status_message("Successfully created $new_name, id: " . $to_model->id);
            }
            else {
                die $self->error_message("Unable to create $new_name model!");
            }
            push @new_models, $to_model;
    }
    $to->assign_models(@new_models); 
    $self->status_message("Monitor this group at: " . $ENV{GENOME_SYS_SERVICES_WEB_VIEW_URL} . "/genome/model-group/status.html?id=" . $to->id);

    return 1;
}

1;
