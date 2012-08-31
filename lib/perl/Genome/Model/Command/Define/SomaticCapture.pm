# FIXME ebelter
#  Long: remove this and all define modeuls to have just one that can handle model inputs
package Genome::Model::Command::Define::SomaticCapture;

use strict;
use warnings;

use Genome;

class Genome::Model::Command::Define::SomaticCapture {
    is => 'Genome::Model::Command::Define::HelperDeprecated',
    has => [
        tumor_model => { 
            is => 'Genome::Model',
            id_by => 'tumor_model_id', 
            doc => 'The tumor model being analyzed'
        },
        tumor_model_id => {
            is => 'Integer',
            is_input => 1,
            doc => 'The tumor model id being analyzed'
        },
        normal_model => { 
            is => 'Genome::Model', 
            id_by => 'normal_model_id', 
            doc => 'The normal model being analyzed'
        },
        normal_model_id => {
            is => 'Integer',
            is_input => 1,
            doc => 'The normal model id being analyzed'
        },
    ],
    has_optional => [
        subject_type => {
            is => 'Text',
            len => 255,
            doc => 'The type of subject all the reads originate from',
            default => 'sample_group',
        },
   ],
};

sub help_synopsis {
    return <<"EOS"
genome model define 
  --subject_name ovc2
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

#    $self->SUPER::execute(@_) or return;

    unless(defined $self->normal_model) {
        $self->error_message("Could not get a model for normal model id: " . $self->normal_model_id);
        return;
    }
    unless(defined $self->tumor_model) {
        $self->error_message("Could not get a model for tumor model id: " . $self->tumor_model_id);
        return;
    }

    #Set up the "subject" of the model
    my $tumor_subject = $self->tumor_model->subject;
    my $normal_subject = $self->normal_model->subject;

    if(($tumor_subject->can('source') || $tumor_subject->can('sample')) and ($normal_subject->can('source') || $normal_subject->can('sample'))) {
        
        my $tumor_source;
        if($tumor_subject->can('source')) {
            $tumor_source = $tumor_subject->source; 
        } else {
            $tumor_source = $tumor_subject->sample->source;
        }
        
        my $normal_source;
        if($normal_subject->can('source')) {
            $normal_source = $normal_subject->source; 
        } else {
            $normal_source = $normal_subject->sample->source;
        }
        
        if($tumor_source eq $normal_source) {
            my $subject = $tumor_source;
            
            #Set up other parameters for call to parent execute()
            $self->subject_id($subject->id);
            $self->subject_class_name($subject->class);
            $self->subject_name($subject->common_name || $subject->name);
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
