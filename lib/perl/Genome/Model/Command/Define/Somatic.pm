# FIXME ebelter
#  Long: remove this and all define modeuls to have just one that can handle model inputs
package Genome::Model::Command::Define::Somatic;

use strict;
use warnings;

use Genome;

class Genome::Model::Command::Define::Somatic {
    is => 'Genome::Model::Command::Define::HelperDeprecated',
    has => [
        tumor_model => { 
            is => 'Genome::Model',
            id_by => 'tumor_model_id', 
            doc => 'The tumor model id being analyzed',
            is_input => 1,
        },
        tumor_model_id => {
            is => 'Integer',
            is_input => 1,
        },
        normal_model => { 
            is => 'Genome::Model', 
            id_by => 'normal_model_id', 
            doc => 'The normal model id being analyzed',
            is_input => 1,
        },
        normal_model_id => {
            is => 'Integer',
            is_input => 1,
        },
        subject_name => {
            is => 'Text',
            is_input => 1,
            is_optional => 1,
            doc => 'Subject name is derived from normal and tumor models and is not necessary as input to somatic models',
        },
    ],
};

sub help_synopsis {
    return <<"EOS"
genome model define 
  --tumor-id 12345
  --normal-id 54321
  --data-directory /gscmnt/somedisk/somedir/model_dir
EOS
}

sub help_detail {
    return <<"EOS"
This defines a new genome model representing the somatic analysis between a normal and tumor model.
EOS
}

sub create {
    my $class = shift;

    my $self = $class->SUPER::create(@_)
        or return;

    return $self;
}

sub type_specific_parameters_for_create {
    my $self = shift;

    my @params = ();

    push @params,
        tumor_model => $self->tumor_model,
        normal_model => $self->normal_model;

    return @params;
}

sub execute {
    my $self = shift;

    unless(defined $self->normal_model) {
        $self->error_message("Could not get a model for normal model id: " . $self->normal_model_id);
        return;
    }
    unless(defined $self->tumor_model) {
        $self->error_message("Could not get a model for tumor model id: " . $self->tumor_model_id);
        return;
    }

    my $tumor_subject = $self->tumor_model->subject;
    my $normal_subject = $self->normal_model->subject;

    if($tumor_subject->can('source') and $normal_subject->can('source')) {
        my $tumor_source = $tumor_subject->source;
        my $normal_source = $normal_subject->source;
        
        if($tumor_source eq $normal_source) {
            my $subject = $tumor_source;
            
            #Set up other parameters for call to parent execute()
            $self->subject($subject);
        } else {
            $self->error_message('Tumor and normal samples are not from same source!');
            return;
        }
    } else {
        $self->error_message('Unexpected subject for tumor or normal model!');
        return;
    }

    # run Genome::Model::Command::Define execute
    my $super = $self->super_can('_execute_body');
    return $super->($self,@_);
}

1;
