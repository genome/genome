package Genome::Model;

use strict;
use warnings;

use Genome;
use Carp;

class Genome::Model {
    is => ['Genome::Notable','Genome::Searchable'],
    is_abstract => 1,
    subclassify_by => 'subclass_name',
    subclass_description_preprocessor => __PACKAGE__ . '::_preprocess_subclass_description',
    id_by => [
        genome_model_id => { is => 'Text', },
    ],
    attributes_have => [
        is_input    => { is => 'Boolean', is_optional => 1, },
        is_param    => { is => 'Boolean', is_optional => 1, },
        is_output   => { is => 'Boolean', is_optional => 1, },
        _profile_default_value => { is => 'Text', is_optional => 1, },
    ],
    has => [
        name => { is => 'Text' },
        type_name => { is => 'Text', via => 'processing_profile' },
        subclass_name => {
            is => 'Text',is_mutable => 0, column_name => 'SUBCLASS_NAME',
            calculate_from => 'processing_profile_id',
            calculate => sub {
                my $pp_id = shift;
                return unless $pp_id;
                my $pp = Genome::ProcessingProfile->get($pp_id);
                unless ($pp) {
                    Carp::croak "Can't find processing profile with ID $pp_id while resolving subclass for model";
                }
                return __PACKAGE__ . '::' . Genome::Utility::Text::string_to_camel_case($pp->type_name);
            },
        },
        _last_complete_build_id => {
            is => 'Number',
            column_name => 'last_complete_build_id',
            doc => 'The last complete build id',
        },
        subject => {
            is => 'Genome::Subject',
            id_by => 'subject_id',
        },
        # FIXME This can be removed once the subject_class_name column is dropped.
        subject_class_name => {
            is => 'Text',
            is_optional => 1,
        },
        processing_profile => {
            is => 'Genome::ProcessingProfile',
            id_by => 'processing_profile_id',
            doc => 'the collection of parameters to be used during the build process'
        },
    ],
    has_optional => [
        limit_inputs_id => {
            is => 'Text',
            column_name => 'LIMIT_INPUTS_TO_ID',
        },
        limit_inputs_rule => {
            is => 'UR::BoolExpr',
            id_by => 'limit_inputs_id',
        },
        user_name => { is => 'Text' },
        creation_date  => { is => 'Timestamp' },
        build_requested => { is => 'Boolean'},
        keep_n_most_recent_builds => { 
            via => 'attributes', to => 'value', is_mutable => 1, 
            where => [ property_name => 'keep_n_most_recent_builds', entity_class_name => 'Genome::Model' ] 
        },
        _last_complete_build_id => { 
            is => 'Number', 
            column_name => 'LAST_COMPLETE_BUILD_ID', 
            doc => 'The last complete build id' ,
        },
        apipe_cron_status => {
            via => 'notes',
            to => 'body_text',
            where => [ header_text => 'apipe_cron_status' ],
            is_mutable => 0,
        },
    ],
    has_optional_many => [
        builds  => {
            is => 'Genome::Model::Build',
            reverse_as => 'model',
            doc => 'Versions of a model over time, with varying quantities of evidence'
        },
        inputs => {
            is => 'Genome::Model::Input',
            reverse_as => 'model',
            doc => 'links to data currently assigned to the model for processing'
        },
        instrument_data => {
            is => 'Genome::InstrumentData',
            via => 'inputs',
            to => 'value',
            is_mutable => 1,
            where => [ name => 'instrument_data' ],
            doc => 'Instrument data currently assigned to the model.'
        },
        projects => {
            is => 'Genome::Project',
            via => 'project_parts',
            to => 'project',
            is_many => 1,
            is_mutable => 1,
            doc => 'Projects that include this model',
        },
        project_parts => {
            is => 'Genome::ProjectPart',
            reverse_as => 'entity',
            is_many => 1,
            is_mutable => 1,
        },
        model_groups => {
            is => 'Genome::ModelGroup',
            via => 'model_bridges',
            to => 'model_group',
            is_mutable => 1
        },
        model_bridges => {
            is => 'Genome::ModelGroupBridge',
            reverse_as => 'model'
        },
    ],
    has_optional_deprecated => [
        auto_assign_inst_data => {
            is => 'Boolean',
        },
        auto_build_alignments => {
            is => 'Boolean',
        },
    ],
    schema_name => 'GMSchema',
    data_source => 'Genome::DataSource::GMSchema',
    table_name => 'GENOME_MODEL',
    doc => 'a versioned data model describing one the sequence and features of a genome'
};

sub define_by {
  # this determines the base class for auto-generated "genome model define XXXX" commands
  # various base classes are available which make different presumptions about presence of instrument-data, etc.
  'Genome::Model::Command::Define::Helper' 
}

# override to do additonal error checking for new profiles
sub __profile_errors__ {
    # This is called whenever a processing profile is created
    my $class = shift;
    my $pp = shift;
    return;
}

# Override in subclasses to have additional stuff appended to the model's default name
# FIXME This will probably go away when the default_name method is overhauled
sub _additional_parts_for_default_name { return; }

# Override in subclasses. Given a list of model inputs missing from a build and a list of build
# inputs missing from the model, should return true if those differences are okay and false otherwise
sub _input_differences_are_ok {
    my $self = shift;
    my @inputs_not_found = @{shift()};
    my @build_inputs_not_found = @{shift()};

    return; #by default all differences are not ok
}

# Override in subclasses. Compares the number of inputs between a build and the model.
sub _input_counts_are_ok {
    my $self = shift;
    my $input_count = shift;
    my $build_input_count = shift;

    return ($input_count == $build_input_count);
}

# Override in subclasses for custom behavior. Updates the model as necessary prior to starting a
# build. Useful for ensuring that the build is incorporating all of the latest information.
# TODO Make sure this is necessary, could be removed
sub check_for_updates {
    return 1;
}

# Override in subclasses, should figure out an appropriate subject for the model and return it
sub _resolve_subject {
    return;
}

# Override in subclasses, should figure out an appropriate processing profile for the model and return it
sub _resolve_processing_profile {
    return;
}

# Override in subclasses
sub _resolve_disk_group_name_for_build {
    # This gets called during the build start process when attempting to create a disk allocation.
    my ($build) = @_;
    return 'info_genome_models';
}

# Override in subclasses
sub _resolve_workflow_for_build {
    # Create a one-step workflow if '_execute_build' is defined.
    my ($self, $build, $optional_lsf_queue) = @_;

    if ($self->can('_execute_build') or $self->processing_profile->can('_execute_build')) { #TODO remove pp._execute_builds
        my $operation_type = Workflow::OperationType::Command->get('Genome::Model::Build::ExecuteBuildWrapper');
        my $resource_requirements = $self->_resolve_resource_requirements_for_build($build);
        $operation_type->lsf_resource($resource_requirements);

        my %opts = (
            name => $build->id . ' all stages',
            input_properties => [ 'build_id' ],
            output_properties => [ 'result' ]
        );

        my $logdir = $build->log_directory;
        if ($logdir =~ /^\/gscmnt/) {
            $opts{log_dir} = $logdir;
        }
 
        my $workflow = Workflow::Model->create(%opts);

        my $operation = $workflow->add_operation(
            name => '_execute_build (' . $self->type_name . ')',
            operation_type => $operation_type, 
        );

        $workflow->add_link(
            left_operation => $workflow->get_input_connector,
            left_property => 'build_id',
            right_operation => $operation,
            right_property => 'build_id'
        );

        $workflow->add_link(
            left_operation => $operation,
            left_property => 'result',
            right_operation => $workflow->get_output_connector,
            right_property => 'result'
        );

        my @e = $workflow->validate;
        die @e unless $workflow->is_valid;

        return $workflow;
    }
    elsif ($self->processing_profile->can("_resolve_workflow_for_build")) {
        return $self->processing_profile->_resolve_workflow_for_build($build);
    }

    my $msg = sprintf(
        "Failed to either implement '_execute_build', or override '_resolve_workflow_for_build', in model '%s' for build '%s'\n",
        $self->__display_name__,
        $build->__display_name__,
    );
    Carp::confess($msg);
}

# Override in subclasses to compose processing profile parameters and build inputs
sub map_workflow_inputs {
    # Construct the input to Workflow::Simple::run_workflow_lsf
    my $self = shift;
    if($self->processing_profile->can('map_workflow_inputs')) {
        return $self->processing_profile->map_workflow_inputs(@_);
    }
    else {
        my ($build) = @_;
        return (build_id => $build->id);
    }
}

# Override in subclasses
sub _resolve_resource_requirements_for_build {
    # This is called during '_resolve_workflow_for_build' to specify the lsf resource requirements of the one-step
    # workflow that is generated from '_execute_build'.
    my ($build) = @_;
    return "-R 'select[model!=Opteron250 && type==LINUX64] rusage[tmp=10000:mem=1000]' -M 1000000";
}

sub _initialize_build {
    # TODO remove call to processing_profile method after moving pp._initialize_builds to the respective model subclasses.
    my ($self, $build) = @_;
    my $pp = $self->processing_profile;
    if ($pp->can('_initialize_build')) {
        return $pp->_initialize_build($build);
    }
    return 1;
}

# Set to true in subclasses if downstream triggering of builds should be enabled
sub _should_trigger_downstream_builds {
    return 0;
}

# Default string to be displayed, can be overridden in subclasses
sub __display_name__ {
    my $self = shift;
    return $self->name . ' (' . $self->id . ')';
}

# Create a model and validate processing profile, subject, etc
sub create {
    my $class = shift;

    # If create is being called directly on this class or on an abstract subclass, SUPER::create will
    # figure out the correct concrete subclass (if one exists) and call create on it.
    if ($class eq __PACKAGE__ or $class->__meta__->is_abstract) {
        return $class->SUPER::create(@_);
    }
    
    my $self = $class->SUPER::create(@_);
    unless ($self) {
        return;
    }

    $self->user_name(Genome::Sys->username) unless $self->user_name;
    $self->creation_date(UR::Context->now);

    $self->_validate_processing_profile;
    $self->_validate_subject;
    $self->_validate_name;

    $self->_verify_no_other_models_with_same_name_and_type_exist;

    # If build requested was set as part of model creation, it didn't use the mutator method that's been
    # overridden. Re-set it here so the required actions take place.
    if ($self->build_requested) {
        $self->build_requested($self->build_requested, 'model created with build requested set');
    }

    return $self;
}

# Delete the model and all of its builds/inputs
sub delete {
    my $self = shift;
    $self->debug_message("Deleting model " . $self->__display_name__);

    my @build_directories;

    for my $input ($self->inputs) {
        $self->debug_message("Deleting model input " . $input->__display_name__);
        my $rv = $input->delete;
        unless ($rv) {
            Carp::confess $self->error_message("Could not delete model input " . $input->__display_name__ .
                " prior to deleting model " . $self->__display_name__);
        }
    }

    for my $build ($self->builds) {
        $self->debug_message("Deleting build " . $build->__display_name__);
        my $rv = $build->delete;
        unless ($rv) {
            Carp::confess $self->error_message("Could not delete build " . $build->__display_name__ .
                " prior to deleting model " . $self->__display_name__);
        }
    }

    return $self->SUPER::delete;
}

# Returns a list of models for which this model is an input
sub to_models {
    my $self = shift;
    my @inputs = Genome::Model::Input->get(value => $self);
    return map { $_->model } @inputs;
}

# Returns a list of models that this model uses as inputs
sub from_models {
    my $self = shift;
    my @inputs = grep { $_->value_class_name->isa('Genome::Model') } $self->inputs;
    return map { $_->value } @inputs;
}

# Returns a list of builds (all statuses) sorted from oldest to newest
sub sorted_builds {
    return shift->builds(-order_by => 'date_scheduled');
}

# Returns a list of succeeded builds sorted from oldest to newest
sub succeeded_builds { return $_[0]->completed_builds; }
sub completed_builds {
    return shift->builds(status => 'Succeeded', -order_by => 'date_completed');
}

# Returns the latest build of the model, regardless of status
sub latest_build {
    my $self = shift;
    my $build_event_iterator = Genome::Model::Event->create_iterator(model_id => $self->id, event_type => 'genome model build', -order_by => '-date_scheduled');
    my $event = $build_event_iterator->next;
    return unless $event;
    return $event->build
}

# Returns the latest build of the model that successfully completed
sub last_succeeded_build { return $_[0]->resolve_last_complete_build; }
sub last_complete_build { return $_[0]->resolve_last_complete_build; }
sub resolve_last_complete_build {
    my $self = shift;
    my $build_event_iterator = Genome::Model::Event->create_iterator(model_id => $self->id, event_type => 'genome model build', event_status => 'Succeeded', -order_by => '-date_completed');
    my $event = $build_event_iterator->next;
    return unless $event;
    return $event->build
}

# Returns a list of builds with the specified status sorted from oldest to newest
sub builds_with_status {
    my ($self, $status) = @_;
    return grep {
        $_->status and
        $_->status eq $status
    } $self->sorted_builds;
}

# Overriding build_requested to add a note to the model with information about who requested a build
sub build_requested {
    my ($self, $value, $reason) = @_;
    # Writing the if like this allows someone to do build_requested(undef)
    if (@_ > 1) {
        my ($calling_package, $calling_subroutine) = (caller(1))[0,3];
        my $default_reason = 'no reason given';
        $default_reason .= ' called by ' . $calling_package . '::' . $calling_subroutine if $calling_package;
        $self->add_note(
            header_text => $value ? 'build_requested' : 'build_unrequested',
            body_text => defined $reason ? $reason : $default_reason,
        );
        return $self->__build_requested($value);
    }
    return $self->__build_requested;
}

# Returns the latest non-abandoned build that has inputs that match the current state of the model
sub current_build {
    my $self = shift;
    my $build_iterator = $self->build_iterator(
        'status not like' => 'Abandoned',
        '-order_by' => '-build_id',
    );
    while (my $build = $build_iterator->next) {
        return $build if $build->is_current;
    }
    return;
}
# Just so current_build_id can be "easily" shown in listers.
sub current_build_id { shift->current_build->id }

# Returns true if no non-abandoned build is found that has inputs that match the current state of the model
sub build_needed {
    return not shift->current_build;
}

# Returns the current status of the model with the corresponding build (if available)
sub status_with_build {
    my $self = shift;
    my ($status, $build);
    if ($self->build_requested) {
        $status = 'Build Requested';
    } elsif ($self->build_needed) {
        $status = 'Build Needed';
    } else {
        $build = $self->current_build;
        $status = ($build ? $build->status : undef);
    }
    return ($status, $build);
}

# Returns the current status of the model
sub status {
    my $self = shift;
    my ($status) = $self->status_with_build;
    return $status;
}

# TODO Clean this up
sub copy {
    my ($self, %overrides) = @_;

    # standard properties
    my %params = ( subclass_name => $self->subclass_name );
    $params{name} = delete $overrides{name} if defined $overrides{name};
    my @standard_properties = (qw/ subject processing_profile auto_assign_inst_data auto_build_alignments /);
    for my $name ( @standard_properties ) {
        if ( defined $overrides{$name} ) { # override
            $params{$name} = delete $overrides{$name};
        }
        elsif ( exists $overrides{$name} ) { # rm undef
            delete $overrides{$name};
        }
        else {
            $params{$name} = $self->$name;
        }
    }

    # input properties
    for my $property ( $self->real_input_properties ) {
        my $name = $property->{name};
        if ( defined $overrides{$name} ) { # override
            my $ref = ref $overrides{$name};
            if ( $ref and $ref eq  'ARRAY' and not $property->{is_many} ) {
                $self->error_message('Cannot override singular input with multiple values: '.Data::Dumper::Dumper({$name => $overrides{$name}}));
                return;
            }
            $params{$name} = delete $overrides{$name};
        }
        elsif ( exists $overrides{$name} ) { # rm undef
            delete $overrides{$name};
        }
        else {
            if ( $property->{is_many} ) {
                $params{$name} = [ $self->$name ];
            }
            else {
                if( defined $self->$name ) {
                    $params{$name} = $self->$name;
                }
            }
        }
    }

    # make we covered all overrides
    if ( %overrides ) {
        $self->error_message('Unrecognized overrides sent to model copy: '.Data::Dumper::Dumper(\%overrides));
        return;
    }

    $params{subject_class_name} = $params{subject}->class; # set here in case subject is overridden

    my $copy = eval{ $self->class->create(%params) };
    if ( not $copy ) {
        $self->error_message('Failed to copy model: '.$@);
        return;
    }

    return $copy;
}

#used for copy above and commands like `genome model input update`
sub real_input_properties {
    my $self = shift;

    my $meta = $self->__meta__;
    my @properties;
    for my $input_property ( sort { $a->property_name cmp $b->property_name } grep { $_->{is_input} or ( $_->via and $_->via eq 'inputs' ) } $meta->property_metas ) {
        my $property_name = $input_property->property_name;
        my %property = (
            name => $property_name,
            is_optional => $input_property->is_optional,
            is_many => $input_property->is_many,
            data_type => $input_property->data_type,
        );

        if($input_property->{is_input}) {
            $property{input_name} = $property_name;
        } else {
            my $where = $input_property->where;
            my %where = @$where;
            $property{input_name} = $where{name};
        }

        if ( $input_property->is_many ) {
            $property{add_method} = 'add_'.$input_property->singular_name,
            $property{remove_method} = 'remove_'.$input_property->singular_name,
        }
        push @properties, \%property;
        next if not $property_name =~ s/_id$//;
        my $object_property = $meta->property_meta_for_name($property_name);
        next if not $object_property;
        $property{name} = $object_property->property_name;
        $property{data_type} = $object_property->data_type;
    }

    return @properties;
}

sub params_for_class {
    my $meta = shift->class->__meta__;

    my @param_names = map {
        $_->property_name
    } sort {
        $a->{position_in_module_header} <=> $b->{position_in_module_header}
    } grep {
        defined $_->{is_param} && $_->{is_param}
    } $meta->property_metas;

    return @param_names;
}

# Called when a build of this model succeeds, requests builds for "downstream" models
# (eg, models that have this model as an input)
sub _trigger_downstream_builds {
    my ($self, $build) = @_;
    return 1 unless $self->_should_trigger_downstream_builds;

    # TODO The observer has to go on build events because that's where build status is 
    # currently stored. Once the status is stored directly on the build itself, the
    # subject_class_name below should be changed to Genome::Model::Build, which should
    # simplify any callbacks as well (since they wouldn't have to check if they are
    # the master event of the build prior to doing anything).
    my @observers = UR::Observer->get(
        subject_class_name => 'Genome::Model::Event::Build',
        aspect => 'event_status',
        note => 'build_success',
    );

    # Only perform this default behavior if no observers for build_success exist for this type
    unless (@observers) {
        my @to_models = $self->to_models;
        for my $model (@to_models) {
            $model->build_requested(
                1,
                'build requested due to successful build ' . $build->id .
                    ' of input model ' . $self->__display_name__
            );
        }
    }

    return 1;
}

# Ensures that processing profile is set. If not, an attempt is made to resolve one before exiting
sub _validate_processing_profile {
    my $self = shift;
    unless ($self->processing_profile) {
        my $pp = $self->_resolve_processing_profile;
        if ($pp and $pp->isa('Genome::ProcessingProfile')) {
            $self->processing_profile($pp);
        }
        else {
            $self->delete;
            Carp::confess "Could not resolve processing profile for model";
        }
    }
    return 1;
}

# Ensures that subject is set. If not, an attempt is made to resolve one before exiting
sub _validate_subject {
    my $self = shift;
    unless ($self->subject) {
        my $subject = $self->_resolve_subject;
        if ($subject and $subject->isa('Genome::Subject')) {
            $self->subject($subject);
        }
        else {
            $self->delete;
            Carp::confess "Could not resolve subject for model";
        }
    }
    return 1;
}

# Ensures that a name is set. If not, a default is used if possible.
sub _validate_name {
    my $self = shift;
    unless ($self->name) {
        my $name = $self->default_model_name;
        if ($name) {
            $self->name($name);
        }
        else {
            $self->delete;
            Carp::confess "Could not resolve default name for model!";
        }
    }
    return 1;
}

# TODO This method should return a generic default model name and be overridden in subclasses.
sub default_model_name {
    my ($self, %params) = @_;

    my $auto_increment = delete $params{auto_increment};
    $auto_increment = 1 unless defined $auto_increment;

    my $name_template = ($self->subject->name).'.';
    $name_template .= 'prod-' if (($self->user_name && $self->user_name eq 'apipe-builder') || $params{prod});

    my $type_name = $self->processing_profile->type_name;
    my %short_names = (
        'genotype microarray' => 'microarray',
        'reference alignment' => 'refalign',
        'de novo assembly' => 'denovo',
        'metagenomic composition 16s' => 'mc16s',
    );
    $name_template .= ( exists $short_names{$type_name} )
    ? $short_names{$type_name}
    : join('_', split(/\s+/, $type_name));

    $name_template .= '%s%s';

    my @parts;
    push @parts, 'capture', $params{capture_target} if defined $params{capture_target};
    push @parts, $params{roi} if defined $params{roi};
    my @additional_parts = eval{ $self->_additional_parts_for_default_name(%params); };
    if ( $@ ) {
        $self->error_message("Failed to get addtional default name parts: $@");
        return;
    }
    push @parts, @additional_parts if @additional_parts;
    $name_template .= '.'.join('.', @parts) if @parts;

    my $name = sprintf($name_template, '', '');
    my $cnt = 0;
    while ( $auto_increment && scalar @{[Genome::Model->get(name => $name)]} ) {
        $name = sprintf($name_template, '-', ++$cnt);
    }

    return $name;
}

# Ensures there are no other models of the same class that have the same name. If any are found, information
# about them is printed to the screen the create fails.
sub _verify_no_other_models_with_same_name_and_type_exist {
    my $self = shift;
    my @models = Genome::Model->get(
        'id ne' => $self->id,
        name => $self->name,
        subclass_name => $self->subclass_name,
    );

    if (@models) {
        my $message = "\n";
        for my $model ( @models ) {
            $message .= sprintf(
                "Name: %s\nSubject Name: %s\nId: %s\nProcessing Profile Id: %s\nSubclass: %s\n\n",
                $model->name,
                $model->subject->name,
                $model->id,
                $model->processing_profile_id,
                $model->subclass_name,
            );
        }
        $message .= sprintf(
            'Found the above %s with the same name and type name.  Please select a new name.',
            Lingua::EN::Inflect::PL('model', scalar(@models)),
        );

        $self->delete;
        Carp::croak $message;
    }

    return 1
}

sub _preprocess_subclass_description {
    my ($class, $desc) = @_;
    my @names = keys %{ $desc->{has} };
    for my $prop_name (@names) {
        my $prop_desc = $desc->{has}{$prop_name};
        # skip old things for which the developer has explicitly set-up indirection
        next if $prop_desc->{id_by};
        next if $prop_desc->{via};
        next if $prop_desc->{reverse_as};
        #next if $prop_desc->{implied_by};

        if ($prop_desc->{is_param} and $prop_desc->{is_input}) {
            die "class $class has is_param and is_input on the same property! $prop_name";
        }

        if (exists $prop_desc->{'is_param'} and $prop_desc->{'is_param'}) {
            $prop_desc->{'via'} = 'processing_profile',
            $prop_desc->{'to'} = $prop_name;
            $prop_desc->{'is_mutable'} = 0;
            $prop_desc->{'is_delegated'} = 1;
            if (defined $prop_desc->{'default_value'}) {
                $prop_desc->{'_profile_default_value'} = delete $prop_desc->{'default_value'};
            }
        }

        if (exists $prop_desc->{'is_input'} and $prop_desc->{'is_input'}) {

            my $assoc = $prop_name . '_association' . ($prop_desc->{is_many} ? 's' : '');
            next if $desc->{has}{$assoc};

            my @where_class;
            if (exists $prop_desc->{'data_type'} and $prop_desc->{'data_type'}) {
                my $prop_class = UR::Object::Property->_convert_data_type_for_source_class_to_final_class(
                    $prop_desc->{'data_type'},
                    $class
                );

                if($prop_class->isa('UR::Value') and !$prop_class->isa('Genome::File::Base')) {
                    push @where_class,
                        value_class_name => $prop_class;
                }
           }

            $desc->{has}{$assoc} = {
                property_name => $assoc,
                implied_by => $prop_name,
                is => 'Genome::Model::Input',
                reverse_as => 'model',
                where => [ name => $prop_name, @where_class ],
                is_mutable => $prop_desc->{is_mutable},
                is_optional => $prop_desc->{is_optional},
                is_many => 1, #$prop_desc->{is_many},
            };

            %$prop_desc = (%$prop_desc,
                via => $assoc,
                to => $class->_resolve_to_for_prop_desc($prop_desc),
            );
        }
    }

    my ($ext) = ($desc->{class_name} =~ /Genome::Model::(.*)/);
    return $desc unless $ext;
    my $pp_subclass_name = 'Genome::ProcessingProfile::' . $ext;

    unless ($desc->{has}{processing_profile}) {
        my $pp_data = $desc->{has}{processing_profile} = {};
        $pp_data->{data_type} = $pp_subclass_name;
        $pp_data->{id_by} = ['processing_profile_id'];

        $pp_data = $desc->{has}{processing_profile_id} = {};
        $pp_data->{data_type} ||= 'Text';
    }

    return $desc;
}

sub _resolve_to_for_prop_desc {
    # TODO This logic was borrowed from the SoftwareResult.pm's _expand_param_input_properties so 
    # when that is refactored, this should also be updated.
    my ($class, $prop_desc) = @_;

    if (exists $prop_desc->{'data_type'} and $prop_desc->{'data_type'}) {
        my $prop_class = UR::Object::Property->_convert_data_type_for_source_class_to_final_class(
            $prop_desc->{'data_type'},
            $class
        );
        if ($prop_class->isa("UR::Value") and ! $prop_class->isa('Genome::File::Base')) {
            return 'value_id';
        } else {
            return 'value';
        }
    } 
    else {
        return 'value_id';
    }
}


my $depth = 0;
sub __extend_namespace__ {
    my ($self,$ext) = @_;

    my $meta = $self->SUPER::__extend_namespace__($ext);
    if ($meta) {
        return $meta;
    }

    $depth++;
    if ($depth>1) {
        $depth--;
        return;
    }

    # If the command class for the model sub type cannot be found, this will create it
    if ( $ext eq 'Command' ) {
        my $create_command_tree = $self->_create_command_tree;
        Carp::confess('Failed to create command tree for '.$self->class.'!') if not $create_command_tree;
    }

    # make a model subclass if the processing profile exists
    # this is deprecated: instead we go the other way and infer the profile from the model
    my $pp_subclass_name = 'Genome::ProcessingProfile::' . $ext;
    my $pp_subclass_meta = eval { $pp_subclass_name->__meta__ };
    if ($pp_subclass_meta and $pp_subclass_name->isa('Genome::ProcessingProfile')) {
        my @pp_delegated_properties = map {
            $_ => { via => 'processing_profile' }
        } $pp_subclass_name->params_for_class;

        my $model_subclass_name = 'Genome::Model::' . $ext;
        my $model_subclass_meta = UR::Object::Type->define(
            class_name => $model_subclass_name,
            is => 'Genome::ModelDeprecated',
            has => \@pp_delegated_properties
        );
        die "Error defining $model_subclass_name for $pp_subclass_name!" unless $model_subclass_meta;
        $depth--;
        return $model_subclass_meta;
    }
    $depth--;
    return;
}

sub _create_command_tree {
    my $self = shift;

    return 1 if $self->__meta__->is_abstract;
    my $command_class_name = $self->class.'::Command';
    my $command_class = eval{ $command_class_name->class; };
    return 1 if $command_class;

    $self->class =~ /^Genome::Model::(\w[\w\d]+)(::)?/;
    my $type_name = $1;
    UR::Object::Type->define(
        class_name => $command_class_name,
        is => 'Command::Tree',
        doc => 'operate on '.Genome::Utility::Text::camel_case_to_string($type_name).' models/builds',
    );

    no strict;
    *{$command_class_name.'::sub_command_category'} = sub { return 'type specific' };
    $command_class = eval{ $command_class_name->class; };

    return $command_class;
}

# diff methods
#
# The build diff methods delegate here so that a pipeline definition doesn't require a hand-written build subclass.
# Only the methods actually used from Build.pm have been migrated.

sub files_ignored_by_build_diff {
    ()
}

1;
