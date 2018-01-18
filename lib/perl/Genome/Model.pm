package Genome::Model;

use strict;
use warnings;

use Genome;
use Genome::Sys::LockProxy qw();
use Carp;
use Scope::Guard;
use UR::Context::Transaction qw(TRANSACTION_STATE_OPEN);

class Genome::Model {
    roles => ['Genome::Role::Notable','Genome::Role::Searchable'],
    table_name => 'model.model',
    is_abstract => 1,
    attributes_have => [
        is_param => {
            is => 'Boolean',
            is_optional => 1,
            doc => 'indicates a value which is part of the processing profile, constant for the model',
        },
        is_input => {
            is => 'Boolean',
            is_optional => 1,
            doc => 'indicates a value which is an input data set for the model: constant per build',
        },
        is_output => {
            is => 'Boolean',
            is_optional => 1,
            doc => 'indicates a value which is produced during the build process',
        },
        is_metric => {
            is => 'Boolean',
            is_optional => 1,
            doc => 'indicates a value which is an interpretation of outputs on the build (can be added post-build)',
        },
        _profile_default_value => {
            is => 'Text',
            is_optional => 1,
            doc => 'on is_param attribute, the default value is stored here, since it is used when making profiles, not making models',
        },
    ],
    subclass_description_preprocessor => 'Genome::Model::_preprocess_subclass_description',
    subclassify_by => 'subclass_name',
    id_by => [
        genome_model_id => {
            # TODO: change to just "id"
            is => 'Text',
            len => 32,
            doc => 'the unique immutable system identifier for a model',
        },
    ],
    has => [
        name => {
            is => 'Text',
            doc => 'the name of the model (changeable)',
        },
        subject => {
            is => 'Genome::Subject',
            id_by => 'subject_id',
            doc => 'the sample/individual/cohort which is being modelled',
        },
        processing_profile => {
            is => 'Genome::ProcessingProfile',
            id_by => 'processing_profile_id',
            doc => 'the collection of parameters to be used during the build process',
        },
        subclass_name => {
            is => 'Text',
            column_name => 'SUBCLASS_NAME',
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

            doc => 'the software subclass for the model in question',
        },
        type_name => {
            is => 'Text',
            via => 'processing_profile',
            doc => 'the name of the type of model (pipeline name)',
        },
        created_by => {
            is => 'Text',
            len => 64,
            doc => 'entity that created the model',
        },
        run_as => {
            is => 'Text',
            len => 64,
            doc => 'username to run builds as',
        },
        creation_date => {
            is => 'DateTime',
            doc => 'the time at which the model was defined',
        },
# These were added by 'ur update classes-from-db' on 22 Sep 2014.  They look like
# old columns that should be removed
#        _user_name => {
#            is => 'Text',
#            len => 64,
#            column_name => 'user_name',
#            is_optional => 1,
#        },
#        build_granularity_unit => {
#            is => 'Text',
#            len => 20,
#            is_optional => 1,
#        },
#        build_granularity_value => {
#            is => 'Integer',
#            len => 4,
#            is_optional => 1,
#        },
#        comparable_normal_model_id => {
#            is => 'Text',
#            len => 32,
#            is_optional => 1,
#        },
#        current_running_build_id => {
#            is => 'Text',
#            len => 32,
#            is_optional => 1,
#        },
#        data_directory => { is => 'Text', len => 1000, is_optional => 1 },
#        is_default => { is => 'Boolean', len => 1, is_optional => 1 },
#        keep => { is => 'Boolean', len => 1, is_optional => 1 },
#        sample_name => { is => 'Text', len => 255, is_optional => 1 },
    ],
    has_optional => [
        user_name => {
            is => 'Text',
            via => '__self__',
            to => 'created_by',
            is_deprecated => 1,
            doc => 'use created_by (or maybe run_as)',
        },
        auto_build => {
            is => 'Boolean',
            doc => 'build automatically when input models rebuild',
            # this is similar to auto_build_alignments, though that flag works on 
            # new instrument data instead of input models
        },
        _build_requested => {
            is => 'Boolean',
            column_name => 'build_requested',
            doc => 'accesses/modifies the "build_requested" value without locking',
        },
        build_requested => {
            via => '__self__',
            to => '_build_requested',
            doc => 'when set to true the system will queue the model for building ASAP',
        },
        _last_complete_build_id => {
            # TODO: change the method with this name to use this property since it is faster
            # nnutter: I disagree, the column should be removed
            # because it was too often out of sync. The method is fairly fast
            # now that it uses an iterator instead of fetching all builds.
            # ssmith: I agree it is more work to write code to get this set correctly.
            # The only issue with that is that having to figure this dynamically is slow.
            # If we do it once when it changes, instead of every time someone wants to check, it could be more DRY.
            is => 'Text',
            column_name => 'LAST_COMPLETE_BUILD_ID',
            doc => 'the last complete build id',
        },
        apipe_cron_status => {
            # This is set in the "genome model build start" command.
            # It is odd for it only to be in the command and not the method, and for it to have a name with apipe in it.
            via => 'notes',
            to => 'body_text',
            is_many => 0,
            where => [ header_text => 'apipe_cron_status' ],
        },
        analysis_project_bridges => {
            is => 'Genome::Config::AnalysisProject::ModelBridge',
            reverse_as => 'model',
            is_many => 1,
        },
        analysis_projects => {
            via => '__self__',
            to => 'analysis_project',
            is_deprecated => 1,
            doc => 'use analysis_project instead',
        },
        analysis_project => {
            is => 'Genome::Config::AnalysisProject',
            via => 'analysis_project_bridges',
            to => 'analysis_project',
            is_many => 0,
        },
        config_profile_item => {
            is => 'Genome::Config::Profile::Item',
            via => 'analysis_project_bridges',
            to => 'config_profile_item',
            is_many => 0,
        },
    ],
    has_many_optional_deprecated => [
        instrument_data => {
            # TODO: the few model types which use instruent data directly should just have:
            #  instrument_data => { is => 'Genome::InstruentData", is_input => 1, is_many => 1 },
            is => 'Genome::InstrumentData',
            via => 'inputs',
            to => 'value',
            is_mutable => 1,
            where => [ name => 'instrument_data' ],
            doc => 'Instrument data currently assigned to the model.',
        },
        model_groups => {
            # TODO: redundant with projects
            is => 'Genome::ModelGroup',
            via => 'model_bridges',
            to => 'model_group',
            is_mutable => 1,
        },
        model_bridges => {
            # TODO: redundant with project_parts
            is => 'Genome::ModelGroupBridge',
            reverse_as => 'model',
        },
    ],
    has_many_optional => [
        builds => {
            is => 'Genome::Model::Build',
            reverse_as => 'model',
            where => [ -order_by => '-created_at' ],
            doc => 'Versions of a model over time, with varying quantities of evidence',
        },
        inputs => {
            is => 'Genome::Model::Input',
            reverse_as => 'model',
            doc => 'links to data currently assigned to the model for processing',
        },
        input_values => {
            via => 'inputs',
            to => 'value',
            doc => 'data currently assigned to the model for processing',
        },

        # TODO: For these to work efficiently we would need last_complete_build to not suck.
        #output_associations => {
        #    via => 'last_complete_build',
        #    to => 'result_users', # rename it in that class to result_associations
        #},
        #outputs => {
        #    via => 'output_associations',
        #    to => 'software_result'
        #},

        project_associations => {
            is => 'Genome::ProjectPart',
            reverse_as => 'part',
            is_mutable => 1,
        },
        projects => {
            is => 'Genome::Project',
            via => 'project_associations',
            to => 'project',
            is_mutable => 1,
            doc => 'Projects that include this model',
        },
    ],
    has_optional_deprecated => [
        _subject_class_name => {
            # TODO: This can be removed once the unused but non-nullable subject_class_name
            # column is dropped.  Subjects all now have a common base class for efficiency.
            # The subject_class_name method is on ::ModelDeprecated and is obsolete.
            is => 'Text',
            column_name => 'SUBJECT_CLASS_NAME',
        },
        subject_class_name => {
            # TODO: this is now read-only, but should go away entirely
            via => 'subject',
            to => 'subclass_name',
        },
        limit_inputs_id => {
            # TODO: this was intended to be a boolexpr ID, but was possibly never used
            is => 'Text',
            column_name => 'LIMIT_INPUTS_TO_ID',
        },
        limit_inputs_rule => {
            # TODO: this was intended to be a boolexpr, but was possibly never used
            is => 'UR::BoolExpr',
            id_by => 'limit_inputs_id',
        },
        auto_assign_inst_data => {
            # this must be here instead of in ::ModelDeprecated becuase it has a db column
            # jeldred noted we set this to false when an apipe-builder model is replaced with a new one
            is => 'Boolean',
        },
        auto_build_alignments => {
            # this must be here instead of in ::ModelDeprecated becuase it has a db column
            # this really means "auto build when instrument data is added"
            is => 'Boolean',
        },
        project_parts => {
            # TODO: the new naming convention was project_associations to prevent confusion
            # between the bridge and the referenced value
            is => 'Genome::ProjectPart',
            reverse_as => 'entity',
            is_mutable => 1,
            is_many => 1,
        },
        _id => {
            # this is the accessor for the column which should become the new primary key
            column_name => 'ID',
        }
    ],
    schema_name => 'GMSchema',
    data_source => 'Genome::DataSource::GMSchema',
    id_generator => '-uuid',
    doc => 'a data model describing the sequence and/or features of a genome',
};

sub define_by {
  # this determines the base class for auto-generated "genome model define XXXX" commands
  # various base classes are available which make different presumptions about presence of instrument-data, etc.
  # TODO: switch to BaseMinimal after making all old models which need ::Helper override this method,
  # including those which are completely auto-generated.
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

# Override in subclasses. Should figure out an appropriate subject for the model and return it.
# TODO: relocate logic in the model define comands from there to here (in subclasses)
sub _resolve_subject {
    return;
}

# Override in subclasses
sub _resolve_disk_group_name_for_build {
    # This gets called during the build start process when attempting to create a disk allocation.
    # TODO: move into the config to have two disk groups, one for "built" things and one for "imported" things
    Genome::Config::get('disk_group_models');
}

# Override in subclasses
sub _resolve_workflow_for_build {
    my ($self, $build) = @_;

    if ($self->can('_execute_build')) {
        # Create a one-step workflow if '_execute_build' is defined.
        my $workflow = Genome::WorkflowBuilder::DAG->create(
            name => $build->workflow_name,
        );

        my $operation = Genome::WorkflowBuilder::Command->create(
            name => '_execute_build (' . $self->type_name . ')',
            command => 'Genome::Model::Build::ExecuteBuildWrapper',
        );
        my $resource_requirements = $self->_resolve_resource_requirements_for_build($build);
        $operation->lsf_resource($resource_requirements);

        $workflow->add_operation($operation);
        $workflow->connect_input(
            input_property => 'build_id',
            destination => $operation,
            destination_property => 'build_id'
        );

        $workflow->connect_output(
            source => $operation,
            source_property => 'result',
            output_property => 'result'
        );

        my $log_dir = $build->log_directory;
        $workflow->recursively_set_log_dir($log_dir);

        return $workflow;
    }
    elsif ($self->processing_profile->can("_resolve_workflow_for_build")) {
        # legacy processing profiles had this logic there
        # TODO: move it to the model subclass instead and remove the PP module
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
    return "-R 'rusage[tmp=10000:mem=1000]' -M 1000000";
}

# Override in subclasses for actions that should be taken at build creation time
sub _initialize_build {
    my ($self, $build) = @_;

    return 1;
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

    my ($bx,%extra) = $class->define_boolexpr(@_);
    if (%extra) {
        # Allow subject_class_name to be specified, though it is
        # no longer used during construction and is a read-only
        # indirect property.
        # TODO: eliminate places which pass this in..
        delete $extra{subject_class_name};
        if (%extra) {
            # If there are still more unknown parameters re-call
            # so this will throw the typical exception minus
            # subject_class_name
            $bx = $class->define_boolexpr($bx->params_list, %extra);
        }
    }

    my $tx = UR::Context::Transaction->begin(commit_validator => sub { 1 });
    my $guard = Scope::Guard->new(sub { local $@; $tx->rollback });

    my $self = $class->SUPER::create($bx);
    unless ($self) {
        return;
    }

    unless ($self->subject) {
        my $subject = $self->_resolve_subject;
        if ($subject and $subject->isa('Genome::Subject')) {
            $self->subject($subject);
        }
        else {
            Carp::confess "Could not resolve subject for model";
        }
    }

    unless ($self->processing_profile) {
        my $reason = eval {
            my $pp_id = $bx->value_for('processing_profile_id');
            unless ($pp_id) {
                return 'No processing profile specified!';
            }

            my $pp = Genome::ProcessingProfile->get($pp_id);
            unless ($pp) {
                return "Specified processing profile ($pp_id) does not exist!";
            }

            my $pp_type = Genome::Utility::Text::string_to_camel_case($pp->type_name);
            my $model_type = ($self->class =~ /^Genome::Model::(.*)/)[0];
            if ($pp_type && $model_type ne $pp_type) {
                return "Processing profile type ($pp_type) does not match model type ($model_type)!";
            }

            return "Missing processing profile; unknown error!";
        };

        Carp::confess $reason || $@;
    }

    for my $m (qw(run_as created_by)) {
        $self->$m(Genome::Sys->username) unless $self->$m;
    }

    unless ($self->name) {
        my $name = $self->default_model_name;
        if ($name) {
            $self->name($name);
        }
        else {
            Carp::confess "Could not resolve default name for model!";
        }
    }

    $self->creation_date(UR::Context->now);

    $self->_verify_no_other_models_with_same_name_and_type_exist;

    # If build requested was set as part of model creation, it didn't use the mutator method that's been
    # overridden. Re-set it here so the required actions take place.
    if ($self->build_requested) {
        $self->build_requested($self->build_requested, 'model created with build requested set');
    }

    if ($self->subject) {
        # TODO: get rid of this as soon as we drop the old database column
        $self->_subject_class_name($self->subject->class);
    }

    # TODO: the column behind this should become the new primary key when we are fully in sync
    $self->_id($self->id);

    $tx->commit() && $guard->dismiss();

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

# Returns a list of succeeded builds sorted from oldest to newest
sub succeeded_builds { return $_[0]->completed_builds; }

sub completed_builds {
    return shift->builds(status => 'Succeeded', -order_by => 'date_completed');
}

# Returns the latest build of the model, regardless of status
sub latest_build {
    my $self = shift;
    my $build_iterator = Genome::Model::Build->create_iterator(model_id => $self->id, -order_by => '-created_at');
    my $build = $build_iterator->next;
    return unless $build;
    return $build
}

# Returns the latest build of the model that successfully completed
sub last_succeeded_build { return $_[0]->resolve_last_complete_build; }
sub last_complete_build { return $_[0]->resolve_last_complete_build; }
sub resolve_last_complete_build {
    my $self = shift;
    my $build_iterator = Genome::Model::Build->create_iterator(model_id => $self->id, status => 'Succeeded', -order_by => '-date_completed');
    my $build = $build_iterator->next;
    return unless $build;
    return $build
}

# Overriding build_requested to add a note to the model with information about who requested a build
sub build_requested {
    my ($self, $value, $reason) = @_;
    # Writing the if like this allows someone to do build_requested(undef)
    if (@_ > 1) {
        if ($value) {
            if ($self->status eq 'Disabled') {
                $self->error_message('Cannot request a build for disabled model %s', $self->__display_name__);
                return;
            } elsif (not $self->analysis_project) {
                $self->error_message('Cannot request a build for model %s without an analysis project', $self->__display_name__);
                return;
            }
        }

        $self->_lock();
        my $default_reason = Carp::shortmess('no reason given');
        $self->add_note(
            header_text => $value ? 'build_requested' : 'build_unrequested',
            body_text => defined $reason ? $reason : $default_reason,
        );
        return $self->__build_requested($value);
    }
    return $self->__build_requested;
}

sub _lock {
    my $self = shift;

    return 1 if $self->{_lock};

    my $model_id = $self->id;
    unless ($ENV{UR_DBI_NO_COMMIT}) {
        my $resource = File::Spec->join('build_requested', $model_id);
        my $lock = $self->{_lock} = Genome::Sys::LockProxy->new(
            resource => $resource,
            scope => 'site',
        )->lock(
            max_try => 30,
            block_sleep => 30,
        );

        die("Unable to acquire the lock to request $model_id. Is something already running or did it exit uncleanly?")
            unless $lock;

        my $cleanup = sub { $lock->unlock };
        my $commit_action = Genome::Sys::CommitAction->create(
            on_commit => $cleanup,
            on_rollback => $cleanup,
        );
    }
}

# Returns the latest non-abandoned build that has inputs that match the current state of the model
sub current_build {
    my $self = shift;
    my $build_iterator = $self->build_iterator(
        'status not like' => 'Abandoned',
        '-order_by' => '-created_at',
    );
    while (my $build = $build_iterator->next) {
        return $build if $build->is_current;
    }
    return;
}

# This is no longer needed since the dot syntax is in place.
# It was made so that current_build.id can be easily shown in listers.
sub current_build_id { shift->current_build->id }

# Returns true if no non-abandoned build is found that has inputs that match the current state of the model
sub build_needed {
    return not shift->current_build;
}

# Returns the current status of the model with the corresponding build (if available)
# Is used only by the "genome model status" command, and underlies the status() property.
sub status_with_build {
    my $self = shift;
    my ($status, $build);

    if ($self->is_disabled) {
        $status = 'Disabled';
    } elsif ($self->build_requested) {
        $status = 'Build Requested';
    } elsif ($self->build_needed) {
        $status = 'Build Needed';
    } else {
        $build = $self->current_build;
        $status = ($build ? $build->status : 'Buildless');
    }
    return ($status, $build);
}

# Returns the current status of the model
sub status {
    my $self = shift;
    my ($status) = $self->status_with_build;
    return $status;
}

sub is_disabled {
    my $self = shift;

    return ($self->config_profile_item and $self->config_profile_item->status eq 'disabled');
}

# Copy model to a new model, overridding some properties
# TODO allow for the filter desc for an input to be passed in
# TODO allow for additions to be passed in
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
            $params{$name} = $self->$name if $self->can($name);
        }
    }

    # input properties
    my @input_properties = $self->real_input_properties;
    for my $property ( @input_properties ) {
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

    my $copy = eval{ $self->class->create(%params) };
    if ( not $copy ) {
        $self->error_message('Failed to copy model: '.$@);
        return;
    }

    # copy the inputs filter desc
    for my $property ( @input_properties ) {
        my $name = $property->{name};
        my @inputs = $self->inputs(name => $name);
        next if not @inputs;
        for my $input ( @inputs ) {
            next if not $input->filter_desc;
            my $copy_input = $copy->inputs(
                name => $name,
                value_class_name => $input->value_class_name, 
                value_id => $input->value_id,
            );
            next if not $copy_input;
            $copy_input->filter_desc( $input->filter_desc );
        }
    }

    return $copy;
}

# used for copy above and commands like `genome model input update`
# TODO: get all models to use the is_input flag instead of having code which is "via" inputs
# and then all of this gets much less complicated.
sub real_input_properties {
    my $self = shift;

    my $meta = $self->__meta__;
    my @properties;
    for my $input_property ( sort { $a->property_name cmp $b->property_name } grep { $_->{is_input} or ( $_->via and $_->via eq 'inputs' ) } $meta->property_metas ) {
        my $property_name = $input_property->property_name;
        next if $property_name eq 'input_values';
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

# TODO This method should return a generic default model name and be overridden in subclasses.
sub default_model_name {
    my ($self, %params) = @_;

    my $name = ($self->subject->name).'.';

    $name .= 'prod-' if (Genome::Sys->user_has_role($self->run_as, 'production') || $params{prod});

    my $type_name = $self->processing_profile->type_name;
    my %short_names = (
        'genotype microarray' => 'microarray',
        'reference alignment' => 'refalign',
        'de novo assembly' => 'denovo',
        'metagenomic composition 16s' => 'mc16s',
        'cwl pipeline' => 'cwl',
    );
    $name .= ( exists $short_names{$type_name} )
    ? $short_names{$type_name}
    : join('_', split(/\s+/, $type_name));

    my @parts = eval{ $self->_additional_parts_for_default_name(%params); };
    if ( $@ ) {
        $self->error_message("Failed to get addtional default name parts: $@");
        return;
    }

    my $suffix = '';
    $suffix .= '.'.join('.', @parts) if @parts;

    return $self->_get_incremented_name($name, $suffix);
}

sub _check_for_existing_model_name {
    my $self = shift;
    my $model_name = shift;
    return scalar @{[Genome::Model->get(name => $model_name)]};
}

sub _get_incremented_name {
    my $self = shift;
    my $model_name = shift;
    my $suffix = shift;
    my $counter = 1;

    my $plain_model_name = $model_name . $suffix;
    unless($self->_check_for_existing_model_name($plain_model_name)) {
        return $plain_model_name;
    }

    my $format_name = sub {
        return sprintf("%s-%s%s", shift, shift, shift);
    };
    while ($self->_check_for_existing_model_name($format_name->($model_name, $counter, $suffix))) {
            $counter++;
    }

    return $format_name->($model_name, $counter, $suffix);
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
            my $prop_class;
            if (exists $prop_desc->{'data_type'} and $prop_desc->{'data_type'}) {
                $prop_class = UR::Object::Property->_convert_data_type_for_source_class_to_final_class(
                    $prop_desc->{'data_type'},
                    $class
                );

                push @where_class,
                    value_class_name => $prop_class;
            }

            my $assoc_property_hash = {
                property_name => $assoc,
                implied_by => $prop_name,
                is => 'Genome::Model::Input',
                reverse_as => 'model',
                where => [ name => $prop_name ],
                is_mutable => $prop_desc->{is_mutable},
                is_optional => $prop_desc->{is_optional},
                is_many => 1, #$prop_desc->{is_many},
            };

            if($prop_class->isa('UR::Value') and !$prop_class->isa('Genome::File::Base')) {
                push @{$assoc_property_hash->{where}}, @where_class;
            }

            $desc->{has}{$assoc} = $assoc_property_hash;

            %$prop_desc = (%$prop_desc,
                via => $assoc,
                to => $class->_resolve_to_for_prop_desc($prop_desc),
            );

            #make a corresponding id property if one doesn't already exist
            if (!$prop_desc->{is_many}) {
                my $id_property = $prop_name . '_id';
                if (!exists($desc->{has}{$id_property})) {
                    $desc->{has}{$id_property} = {
                        property_name => $id_property,
                        implied_by => $prop_name,
                        via => $assoc,
                        to => 'value_id',
                        where => \@where_class,
                        is_mutable => $prop_desc->{is_mutable},
                        is_optional => $prop_desc->{is_optional},
                        is_many => 0,
                    };
                }
            }
        }

    }

    my ($ext) = ($desc->{class_name} =~ /Genome::Model::(.*)/);
    return $desc unless $ext;
    my $pp_subclass_name = 'Genome::ProcessingProfile::' . $ext;

    unless ($desc->{has}{processing_profile}) {
        my $pp_data = $desc->{has}{processing_profile} = {};
        $pp_data->{data_type} = $pp_subclass_name;
        $pp_data->{id_by} = ['processing_profile_id'];
        $pp_data->{doc} = Genome::Model->__meta__->property('processing_profile')->doc;

        my $pp_id_data = $desc->{has}{processing_profile_id};
        unless ($desc->{has}{processing_profile_id}) {
            $pp_id_data = $desc->{has}{processing_profile_id} = {};
            $pp_id_data->{implied_by} = 'processing_profile';
            $pp_id_data->{doc} = Genome::Model->__meta__->property('processing_profile_id')->doc;
        }
        $pp_id_data->{data_type} ||= 'Text';
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
        my $create_command_tree = $self->_extend_namespace_with_command_tree;
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

sub _extend_namespace_with_command_tree {
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

sub files_ignored_by_build_diff { () }

#Does this model type require mapping of subjects to work?
#Used by Analysis Project configuration to "pair" somatic samples or otherwise aggregate them for analysis
sub requires_subject_mapping { return 0; }

sub sudo_wrapper { '/usr/local/bin/sudo-genome-build-start' }

sub should_run_as {
    my $self = shift;
    return unless $self->_can_run_as();
    return (Genome::Sys->username ne $self->run_as);
}

sub _can_run_as {
    my $self = shift;

    my $sudo_wrapper = sudo_wrapper();
    unless ( -x $sudo_wrapper ) {
        return;
    }

    return 1;
}

1;
