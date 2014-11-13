package Genome::ModelGroup::Command::Create;

use strict;
use warnings;

use Genome;

class Genome::ModelGroup::Command::Create {
    is => 'Command::V2',
    has => [
        name => {
            is => 'Text',
            len => 255, 
            shell_args_position => 1,
            doc => 'A name for the model group.', 
        },
        models => {
            is => 'Genome::Model',
            is_many => 1,
            is_optional => 1,
            shell_args_position => 2,
            doc => 'Model(s) to add to group. Resolved from command line via text string.',
        },
    ],
};

sub help_brief {
    return "create a new model-group";
}

sub help_detail {
    return "create a new model-group";
}

sub help_synopsis {
    return <<HELP;
    genome model-group create "Example Group Name" [MODELS]

  Example:
    genome model-group create "TCGA-UCEC-build36-SSv2broad-ckandoth" 2868304369 2868304371 2868304372 ...
HELP
}

sub execute {
    my $self = shift;

    my @models = $self->models;

    my $model_group = Genome::ModelGroup->create(
        name => $self->name,
    );
    
    if ( not $model_group ) {
        $self->error_message('Failed to create model group for name: '.$self->name);
        return;
    }
    $self->status_message('Created model group:');
    $self->status_message('ID: ' . $model_group->id . ', NAME: ' . $model_group->name);
    
    if ( @models ) {
        $self->status_message('Assigning '.scalar(@models).' to group...');
        $model_group->assign_models(@models);
    }
    
    return 1;
}

1;

