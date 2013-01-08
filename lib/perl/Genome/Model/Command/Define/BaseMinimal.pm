package Genome::Model::Command::Define::BaseMinimal;

# the other base classes for model define are:
#   Helper.pm           the latest and greatest ...but expects instrument data to be part of it, generated dynamically by ::Define
#   Base.pm             an alternative to Helper
#   HelperDeprecated.pm extends helper with deprecated methods 

# This is a subset of Helper which does not presume:
# 1. that there are instrument_data on the model
# 2. that the user directly specifies the subject
# 3. that the user directly specifies the processing_profile
#
# It does, however, allow the user to set any inputs like the autogen
# Helper subclasses.
#
# Best practice is to write a subclass of this, and:
# - add instrument_data, processing_profile, and subject if/as needed.
# - add other params to resolve those dynamically.

use strict;
use warnings;

use Genome;
use File::Path;
use Carp 'confess';

class Genome::Model::Command::Define::BaseMinimal {
    is => 'Command::V2',
    is_abstract => 1,
    subclass_description_preprocessor => __PACKAGE__ . '::_preprocess_subclass_description',
    has_optional_input => [
        add_to_projects => {
            is => 'Genome::Project',
            is_many => 1,
            doc => 'add the new model to these projects (model groups)',
        },        
    ],
    has_output => [
        model => {
            is => 'Genome::Model',
            id_by => 'result_model_id',
            is_transient => 1,
            doc => 'a new genome model',
        },
    ],
    doc => "define a new genome model"
};

sub _target_class_name {
    my $class = shift;
    $class = $class->class;
    $class =~ s/::Command::Define::/::/;
    return $class;
}

sub _resolve_subject_from_inputs { return }

sub _preprocess_subclass_description {
    my ($this_class_name, $desc) = @_;
    my $cmd_subclass_name = $desc->{class_name};
    my $model_subclass_name = $cmd_subclass_name;
    $model_subclass_name =~ s/::Command::Define::/::/;

    my $model_subclass_meta = UR::Object::Type->get($model_subclass_name);
    if ($model_subclass_meta and $model_subclass_name->isa('Genome::Model')) {
        # determine the exact base class dynamically (some subclass of this one)
        my @inheritance = ($this_class_name);
        if ($model_subclass_name->can("define_by")) {
            my $base_class = $model_subclass_name->define_by;
            $base_class->class;
            unless ($base_class->isa(__PACKAGE__)) {
                die "$model_subclass_name has a define_by value of $base_class, but this is not a subclass of BaseMinimal??";
            }   
            @inheritance = ($base_class);
        }
            
        # add inputs to the define command where the model has inputs

        # ...but do not support certain inputs at definition time
        my %suppress;
        if ($cmd_subclass_name->can('_suppress_inputs')) {
            %suppress = map { $_ => 1 } $cmd_subclass_name->_suppress_inputs();
        }
        
        my @p = $model_subclass_meta->properties();
        my $has = $desc->{has};
        for my $p (@p) {
            my $name = $p->property_name;
            next if $has->{$name};
            next if $suppress{$name};
            if (($p->can("is_input") and $p->is_input) or $name =~ /^(processing_profile|processing_profile_id|name)$/) {
                my %data = %{ UR::Util::deep_copy($p) };
                for my $key (keys %data) {
                    delete $data{$key} if $key =~ /^_/;
                }
                delete $data{id};
                delete $data{db_committed};
                delete $data{via};
                delete $data{to};
                $data{is_input} = 1;
                if ($name eq 'name') {
                    # support --name and --model-name, giving prefernece to the former
                    if (grep { $_->can("model_name") or $_->can("name") } @inheritance) {
                        # a base class specifies model_name or name
                        next;
                    }
                    if ($has->{model_name} or $has->{name}) {
                        # this desc specifies model_name or name
                        next;
                    }
                    $data{doc} = 'a friendly name for the new model (changeable)';
                    $data{is_optional} = 1;
                }
                $has->{$name} = \%data;
            }
        }

        # add an input for model_name if it is not present

        # add an input for subject unless it is a calculated property (TODO)
    }
    
    return $desc;
}

sub help_brief {
    my $self = shift;
    my $msg;
    my $model_subclass = $self->_target_class_name;
    if ($model_subclass) {
        my $model_type = eval { $model_subclass->__meta__->property('processing_profile')->data_type->_resolve_type_name_for_class };
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
    my $target_class_name = $self->_target_class_name;
    if ($target_class_name->can("_help_synopsis_for_model_define")) {
        return $target_class_name->_help_synopsis_for_model_define;
    }
    elsif ($target_class_name->can("_help_synopsis")) {
        return $target_class_name->_help_synopsis;
    }
    #my $model_type = $target_class_name->__meta__->property('processing_profile')->data_type->_resolve_type_name_for_class;
    #$model_type =~ s/[_ ]/-/g;
    #my @show = "  genome model define $model_type";
    my $command_name = $self->command_name;
    my @show = "  $command_name";
    if ($self->can("model_name")) {
        push @show, "    --model-name test1";
    }
    elsif ($self->can("name")) {
        push @show, "    --name test1";
    }
    if ($self->can("subject")) {
        push @show, "    --subject TEST-patient1-sample1";
    }
    push @show, "    --processing-profile name='my processing profile'";
    my $txt = join(" \\\n",@show) . "\n";
    return $txt;
}

sub help_detail {
    my $self = shift;
    my $target_class_name = $self->_target_class_name; 
    if ($target_class_name->can("_help_detail_for_model_define")) {
        return $target_class_name->_help_detail_for_model_define;
    }
    elsif ($target_class_name->can("_help_detail")) {
        return $target_class_name->_help_detail;
    }
    return <<"EOS"
This defines a new genome model for the specified subject, using the 
specified processing profile.
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
    
    unless ($self->_validate_inputs) {
        confess "Could not validate inputs!";
    }

    unless ($self->subject) {
        my $subject = $self->_resolve_subject_from_inputs;
        if ($subject) {
            $self->subject($subject);
        }
        # if this fails we still try to let the constructor do it...
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
        $self->type_specific_parameters_for_create,
    );

    # something odd is happening when passing in the object...
    if (my $p = delete $params{processing_profile}) {
        $params{processing_profile_id} = $p->id
    }

    if (my @projects = $self->add_to_projects) {
        $params{projects} = \@projects;
    }

    #print Data::Dumper::Dumper(\%params);
    my $target_class = $self->_target_class_name;
    my $model = $target_class->create(%params);
    unless ($model) {
        confess "Could not create a model!";
    }

    unless ($model->subject) {
        $model->delete;
        Carp::confess("Not given subject and could not derive subject from inputs!");
    }

    $self->model($model);
    $self->display_model_information($model);
    return 1;
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

