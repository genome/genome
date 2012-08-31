package Genome::Model::Command::Define::Convergence;

use strict;
use warnings;

use Genome;

class Genome::Model::Command::Define::Convergence {
    is => 'Genome::Model::Command::Define::HelperDeprecated',
    has => [
        model_group_id => {
            is => 'Text',
            doc => 'The id of the model group for which to create a convergence model'
        },
        model_group => {
            is => 'Genome::ModelGroup',
            id_by => 'model_group_id',
        },
    ],
    has_optional => [
        model_name => {
            is => 'Text',
            len => 255,
            doc => 'User meaningful name for this model (defaults to the model group name with "_convergence" appended',
        },
        processing_profile_name => {
            is => 'Text',
            doc => 'identifies the processing profile by name',
            default => 'convergence default',
        },
        auto_build_alignments => { #TODO Yeah, they're not really "alignments", but that's what the parameter is in the parent class
            is => 'Boolean',
            doc => 'If true, new builds will automatically be launched when the underlying model group changes.',
            default_value => 1,
        },
   ],
};

sub help_synopsis {
    return <<"EOS"
genome model define convergence 
  --model-group-id 242 
  --data-directory /gscmnt/somedisk/somedir/model_dir
EOS
}

sub help_detail {
    return <<"EOS"
This defines a new genome model representing the harmonic convergence analysis for a group of models.
EOS
}

sub create {
    my $class = shift;
    
    my $self = $class->SUPER::create(@_)
        or return;

    unless(defined $self->model_group) {
        $self->error_message('No ModelGroup found for id: ' . $self->model_group_id);
        return;
    }

    unless($self->model_name) {
        $self->model_name($self->model_group->name . '_convergence');
    }

    my $subject = $self->model_group->infer_group_subject;
    $self->subject($subject);

    return $self;
}

sub execute {
    my $self = shift;
    
    #execute for commands is renamed after some Command.pm magic
    $self->SUPER::_execute_body(@_) or return;

    # get the model created by the super
    my $model = Genome::Model->get($self->result_model_id);

    unless($model){
        $self->error_message('Could not get model from base define command.');
        return;
    }
    
    my $set_group_cmd = Genome::Model::Command::Input::Update->create(
        models => [$model],
        name => 'group',
        value => $self->model_group->id,
    );
    $set_group_cmd->dump_status_messages(1);
    unless($set_group_cmd->execute) { 
        $self->error_message('Could not set group for model.');
        return;
    }
    
    return 1;
}

1;
