package Genome::Model::Command::SetDefault;

use strict;
use warnings;

use Genome;

class Genome::Model::Command::SetDefault {
    is => ['Command'],
    has => [
        model           => { is => 'Genome::Model', id_by => 'model_id' },
        model_id        => { is => 'Text', doc => 'identifies the genome model by id', shell_args_position => 1 },
    ],
    has_optional => [
        clear           => { is => 'Boolean', default => '0', doc => 'Indicate the provided model should no longer be considered the default.' },
        replace         => { is => 'Boolean', default => '0', doc => 'Set the provided model as default even if it means replacing another model already set as the default.'}
    ]
};

sub help_brief {
    return "A simple tool for setting/clearing default models.";
}

sub help_synopsis {
    my $self = shift;
    return <<"EOS"
    genome model set-default 12345
    genome model set-default --clear 12345
    genome model set-default --replace 12345
EOS
}

sub help_detail {                           
    return <<EOS 
A subject may have many models associated with it.  This tool provides a mechanism by which to control which of these models is considered the "default" for purposes of analysis.
EOS
}

sub execute {
    my $self = shift;
    
    my $model = $self->model;
    
    unless ($model) {
        $self->error_message('Could not locate a model with the given parameters.');
        return;
    }
    
    $self->status_message('Model: ' .  $model->name . ' (' . $model->id . ')');
    
    if($self->clear) {
        #Can skip most of the checks--this can't lead to too many default models
        $model->is_default(undef);
        $self->status_message('Cleared default flag on model.');
        return 1;
    }
    
    unless ($model->subject) {
        $self->error_message('Could not discern a subject for the provided model.');
        return;
    }
    
    $self->status_message( 'Found subject for model: ' . $model->subject_name . ' (' . $model->subject_id . ')');
    
    #Find all other models with the same subject
    my @other_models = Genome::Model->get( subject_id => $model->subject_id, subject_class_name => $model->subject_class_name );
    
    for my $other_model (@other_models) {
        next if $model->id eq $other_model->id;
        
        #We're looking for any existing models that might be the default.
        next unless $other_model->is_default;
        
        if($self->replace) {
            $self->status_message('Clearing default flag on model: ' .  $other_model->name . ' (' . $other_model->id . ')');
            $other_model->is_default(undef);
        } else {
            $self->status_message('Found existing default model for this subject: ' .  $other_model->name . ' (' . $other_model->id . ').  Use --replace if you want to change the default model from that model to the provided one.' );
            return;
        }
    }
    
    $model->is_default(1);
    $self->status_message('Set model as default for subject.');
    
    return 1;
}

1;
