package Genome::Model::Command::Define::Helper;

use strict;
use warnings;

use Genome;
use File::Path;
use Carp 'confess';

class Genome::Model::Command::Define::Helper {
    is => 'Command::V2',
    is_abstract => 1,
    has => [
        processing_profile => {
            is => 'Genome::ProcessingProfile',
            id_by => 'processing_profile_id',
            doc => 'Processing profile to be used by model, can provide either a name or an id',
        },
    ],
    has_optional => [
        subject => {
            is => 'Genome::Subject',
            id_by => 'subject_id',
            doc => 'Subject for the model, can provide either a name or an id. If instrument data is provided and this is not, ' .
                'an attempt will be made to resolve it based on the provided instrument data'
        },
        model_name => {
            is => 'Text',
            doc => 'User meaningful name for this model, a default is used if none is provided',
        },
        auto_assign_inst_data => {
            is => 'Boolean',
            default_value => 0,
            doc => 'Assigning instrument data to the model is performed automatically',
        },
        auto_build_alignments => {
            is => 'Boolean',
            default_value => 1,
            doc => 'The building of the model is performed automatically',
        },
        result_model_id => {
            is => 'Number',
            is_transient => 1,
            doc => 'Stores the ID of the newly created model, useful when running this command from a script',
        },
    ],
    has_many_optional => [
        instrument_data => {
            is => 'Genome::InstrumentData',
            doc => 'Instrument data to be assigned to the model, can provide a query to resolve, a list of ids, etc'
        },
        groups => {
            is => 'Genome::ModelGroup',
            doc => 'Model groups to put the newly created model into',
        },        
        projects => {
            is => 'Genome::Project',
            doc => 'Projects for the model.',
        },
    ],
};

sub help_brief {
    my $self = shift;
    my $msg;
    my $model_subclass = eval {$self->_target_class_name };
    if ($model_subclass) {
        my $model_type = $model_subclass->__meta__->property('processing_profile')->data_type->_resolve_type_name_for_class;
        if ($model_type) {
            $msg = "define a new $model_type genome model";
        }
        else {
            $msg = "define a new genome model";
        }
    }
    else {
        $msg = 'define a new genome model'
    }
    return $msg;
}

sub help_synopsis {
    my $self = shift;
    my $target_class_name = $self->class;
    $target_class_name =~ s/::Command::Define::/::/;
    if ($target_class_name->can("_help_synopsis_for_model_define")) {
        return $target_class_name->_help_synopsis_for_model_define;
    }
    elsif ($target_class_name->can("_help_synopsis")) {
        return $target_class_name->_help_synopsis;
    }
    return <<"EOS"
genome model define reference-alignment 
  --model-name test5
  --subject ley_aml_patient1_tumor
  --processing-profile nature_aml_08
EOS
}

sub help_detail {
    my $self = shift;
    my $target_class_name = $self->class;
    $target_class_name =~ s/::Command::Define::/::/;
    if ($target_class_name->can("_help_detail_for_model_define")) {
        return $target_class_name->_help_detail_for_model_define;
    }
    elsif ($target_class_name->can("_help_detail")) {
        return $target_class_name->_help_detail;
    }
    return <<"EOS"
This defines a new genome model for the specified subject, using the 
specified processing profile.

The first build of a model must be explicitly requested after the model is defined.
EOS
}

# Should be overridden in subclasses
sub type_specific_parameters_for_create {
    my $self = shift;
    my $model_class = $self->class;
    $model_class =~ s/::Command::Define::/::/;
    my @p = 
        map { 
            my $meta = $_;
            my $name = $meta->property_name;
            my @values = grep { defined $_ } $self->$name;
            if (@values == 0) {
                ()
            }
            elsif ($meta->is_many or @values > 1) {
                ($name => \@values);

            }
            else {
                ($name => $values[0]);
            }
        }
        grep { $model_class->can($_->property_name) }
        grep { $_->can("is_input") and $_->is_input }
        $self->__meta__->properties();
    
    return @p 
}

sub execute {
    my $self = shift;

    my $processing_profile = $self->validate_processing_profile;
    unless ($processing_profile) {
        confess "Could not validate processing profile!";
    }

    unless ($self->subject) {
        my $subject = $self->deduce_subject_from_instrument_data;
        if ($subject) {
            $self->subject($subject);
        }
    }

    my %params = (
        (   
            $self->subject ? 
            (
                subject_id => $self->subject->id,
                subject_class_name => $self->subject->class,
            )
            :
            ()
        ),
        processing_profile_id => $self->processing_profile->id,
        name => $self->model_name,
        auto_assign_inst_data => $self->auto_assign_inst_data,
        auto_build_alignments => $self->auto_build_alignments,
        $self->type_specific_parameters_for_create,
    );
    $params{instrument_data} = [$self->instrument_data] if $self->instrument_data;
    $params{model_groups} = [$self->groups] if $self->groups;
    $params{projects} = [$self->projects] if $self->projects;
            
    my $model = Genome::Model->create(%params);
    unless ($model) {
        confess "Could not create a model!";
    }

    unless ($model->subject) {
        my @instrument_data = $self->instrument_data;
        $model->delete;
        if (@instrument_data) {
            confess "Not given subject and could not derive subject from instrument data!";
        }
        else {
            confess "No instrument data provided, cannot deduce subject!";
        }
    }

    $self->result_model_id($model->id);
    $self->display_model_information($model);
    return 1;
}

sub deduce_subject_from_instrument_data {
    my $self = shift;

    my @instrument_data = $self->instrument_data;
    unless (@instrument_data) {
        return; 
    }

    my @samples = map { $_->sample } @instrument_data;
    unless (grep { $_->id ne $samples[0]->id } @samples) {
        return $samples[0];
    }

    my @individuals = map { $_->source } @samples;
    unless (grep { $_->id ne $individuals[0]->id } @individuals) {
        return $individuals[0];
    }

    my @taxons = map { $_->taxon } @individuals;
    unless (grep { $_->id ne $taxons[0]->id } @taxons) {
        return $taxons[0];
    }

    return;
}

sub display_model_information {
    my ($self, $model) = @_;

    $self->status_message("Created model:");
    my $list = Genome::Model::Command::List->create(
        filter => 'id=' . $model->id,
        show => join(',', $self->listed_params),
        style => 'pretty',
    );
    return $list->execute;
}

sub listed_params {
    return qw/ id name subject.name subject.subject_type processing_profile_id processing_profile.name /;
}

sub validate_processing_profile {
    my $self = shift;
    unless ($self->compare_pp_and_model_type) {
        confess 'Model and processing profile types do not match!';
    }
    return 1;
}

# TODO This may not be necessary. Even if it is, there's probably a better way to do it.
sub compare_pp_and_model_type {
    my $self = shift;

    # Determine the subclass of model being defined
    my $model_subclass = $self->class;
    my $package = "Genome::Model::Command::Define::";
    $model_subclass =~ s/$package//;
    
    # Determine the subclass of the processing profile
    my $pp = Genome::ProcessingProfile->get(id => $self->processing_profile->id);
    unless($pp){
        $self->error_message("Couldn't find the processing profile identified by the #: " . $self->processing_profile->id);
        die $self->error_message;
    }
    my $pp_subclass = $pp->subclass_name;
    $pp_subclass =~ s/Genome::ProcessingProfile:://;
    
    unless ($model_subclass eq $pp_subclass) {
        my ($shortest, $longest) = ($model_subclass, $pp_subclass);
        ($shortest, $longest) = ($longest, $shortest) if length $pp_subclass < length $model_subclass;
        unless ($longest =~ /$shortest/) {
            $self->error_message("Model subclass $model_subclass and ProcessingProfile subclass $pp_subclass do not match!");
            confess;
        }
    }

    return 1;
}

1;

