package Genome::Model::Build;

use strict;
use warnings;

use Genome;

use Carp;
use Data::Dumper 'Dumper';
use File::stat;
use File::Path;
use File::Find 'find';
use File::Next;
use File::Basename qw/ dirname fileparse /;
use Regexp::Common;
use Workflow;
use YAML;
use Date::Manip;

use Genome::Utility::Email;

class Genome::Model::Build {
    is => ['Genome::Notable','Genome::Searchable'],
    type_name => 'genome model build',
    table_name => 'model.build',
    is_abstract => 1,
    subclassify_by => 'subclass_name',
    subclass_description_preprocessor => __PACKAGE__ . '::_preprocess_subclass_description',
    id_by => [
        # TODO: change to just "id"
        build_id => { is => 'Text', len => 64 },
    ],
    attributes_have => [
        is_input    => { is => 'Boolean', is_optional => 1, },
        is_param    => { is => 'Boolean', is_optional => 1, },
        is_output   => { is => 'Boolean', is_optional => 1, },
        is_metric   => { is => 'Boolean', is_optional => 1 },
    ],
    has => [
        is_last_complete => {
            is => 'Boolean',
            calculate_from => ['id','model'],
            calculate => q|my $build = $model->last_complete_build; return $build && $build->id eq $id|,
            doc => 'true for any build which is the last complete bulid for a model',
        },
        subclass_name => {
            is => 'Text',
            len => 255,
            is_mutable => 0,
            column_name => 'SUBCLASS_NAME',
            calculate_from => ['model_id'],
            # We subclass via our model's type_name (which is via it's processing profile's type_name)
            calculate => sub {
                my($model_id) = @_;
                return unless $model_id;
                my $model = Genome::Model->get($model_id);
                Carp::croak("Can't find Genome::Model with ID $model_id while resolving subclass for Build") unless $model;
                return __PACKAGE__ . '::' . Genome::Utility::Text::string_to_camel_case($model->type_name);
            }
        },
        data_directory          => { is => 'Text', len => 1000, is_optional => 1 },
        model                   => { is => 'Genome::Model', id_by => 'model_id' },
        type_name               => { via => 'model', is => 'Text' },
        subject                 => { via => 'model', is => 'Genome::Subject' },
        processing_profile      => { via => 'model', is => 'Genome::ProcessingProfile' },
        run_by                  => { via => 'creation_event', to => 'user_name' },
        status                  => { via => 'creation_event', to => 'event_status', is_mutable => 1 },
        date_scheduled          => { via => 'creation_event', to => 'date_scheduled', },
        date_completed          => { via => 'creation_event', to => 'date_completed' },

        # this is the one event type we tentatively intended to keep for builds
        # even fully workflowified, it is possibly still a candidate for retention, logging things like input changes, etc.
        creation_event          => { is => 'Genome::Model::Event', reverse_as => 'build', where => [ event_type => 'genome model build' ], is_many => 1, is_constant => 1},
        _events                 => { is => 'Genome::Model::Event', reverse_as => 'build', is_many => 1 },
    ],
    has_optional => [
        _newest_workflow_instance => {
            is => 'Workflow::Operation::Instance',
            is_calculated => 1,
            calculate => q{ return $self->newest_workflow_instance(); }
        },
        disk_allocation   => { 
            is => 'Genome::Disk::Allocation', 
            calculate_from => [ 'class', 'id' ],
            calculate => q(
                my $disk_allocation = Genome::Disk::Allocation->get(
                                        owner_class_name => $class,
                                        owner_id => $id,
                                    );
                return $disk_allocation;
            ) 
        },
        software_revision => { 
            is => 'Text', 
            len => 1000 
        },
    ],
    has_many_optional => [
        input_values => {
            via => 'input_associations',
            to => 'value',
            doc => "The values associated with a build's input associations.",
        },
        input_associations => {
            is => 'Genome::Model::Build::Input',
            reverse_as => 'build',
            doc => 'Links between a build and its input values, including the specification of which input the value satisfies.'
        },
        downstream_build_associations => {
            is => 'Genome::Model::Build::Input',
            reverse_as => '_build_value',
            doc => 'links to models which use this model as an input', 
        },
        downstream_builds => {
            is => 'Genome::Model::Build',
            via => 'downstream_build_associations',
            to => 'build',
            doc => 'models which use this model as an input',
        },

        instrument_data  => {
            is => 'Genome::InstrumentData',
            via => 'input_associations',
            to => 'value',
            is_mutable => 1,
            where => [ name => 'instrument_data' ],
            doc => 'Instrument data assigned to the model when the build was created.'
        },
        instrument_data_associations => {
            is => 'Genome::Model::Build::Input',
            reverse_as => 'build',
            where => [ name => 'instrument_data' ],
        },
        results => {
            # TODO rename to maybe outputs?
            is => 'Genome::SoftwareResult',
            via => 'result_users',
            to => 'software_result',
        },
        result_associations => {
            # TODO rename to maybe result_associations or output_associations?
            # This implies that we are bridging to something using the build.
            is => 'Genome::SoftwareResult::User',
            reverse_as => 'user',
        },
        attribute_associations => { 
            is => 'Genome::MiscAttribute', 
            reverse_as => '_build', 
            where => [ entity_class_name => 'Genome::Model::Build' ] 
        },
        metrics => { 
            is => 'Genome::Model::Metric', 
            reverse_as => 'build',
            doc => 'Build metrics' 
        },
        projects => {
            # TODO: shoudln't this be Genome::Project?  Genome::Site::TGI::Project is used by the sync cron only, right?
            is => 'Genome::Site::TGI::Project', 
            via => 'model' 
        },
        work_orders => {
            # TODO: is this relationship correct?
            # workorders should be more granular than projects
            is => 'Genome::WorkOrder', 
            via => 'projects' 
        },
    ],
    has_optional_deprecated => [
        # unnecessary now that the dot syntax is supported on the command line, and also in get():
        processing_profile_id   => { via => 'model' },
        processing_profile_name => { via => 'model' },
        subject_id              => { via => 'model' },
        subject_name            => { via => 'subject', to => 'name' },
        model_id                => { is => 'NUMBER', implied_by => 'model', constraint_name => 'GMB_GMM_FK' },
        model_name              => { via => 'model', to => 'name' },
        
        # this should be moved down into the classes which actually use it
        region_of_interest_set_name => {
            is => 'Text',
            is_many => 1,
            is_mutable => 1,
            via => 'inputs',
            to => 'value_id',
            where => [ name => 'region_of_interest_set_name', value_class_name => 'UR::Value' ],
        },
        
        # poorly named 
        the_master_event => { 
            is => 'Genome::Model::Event', 
            via => 'creation_event', 
            to => '__self__' 
        }, 
    ],
    has_many_optional_deprecated => [
        inputs => { is => 'Genome::Model::Build::Input', via => '__self__', to => 'input_associations' },
        
        # the_* ?
        # I don't believe anythign which uses events actually needs these
        the_events => { is => 'Genome::Model::Event', to => '_events', via => '__self__' },
        the_events_statuses => { is => 'Text', to => '_events', via => 'event_status' },
        
        # these are ambiguosuly named and are still present only for backward compatibility
        # these returns the "bridge", but are named for the thing on the other side of the bridge
        # now there is "*_associations", which are clearly a bridge,
        # and "*_values", which are unambiguously the value across the bridge
        attributes => { is => 'Genome::MiscAttribute', via => 'attribute_associations', to => '__self__', },
        instrument_data_inputs => { is => 'Genome::Model::Build::Input', to => 'instrument_data_associations', via => '__self__', },
        result_users => {
            # this is particularly confusing, because it is giving the user
            # of the builds's results, but the cases where the build itself
            # is the user of the result (ie: use $build->result_associations)
            is => 'Genome::SoftwareResult::User',
            to => 'result_associations',
            via => '__self__',
        },
       
        # unnecessary now that the dot syntax is supported on the command line, and also in get():
        instrument_data_ids => { via => 'instrument_data', to => 'id', is_many => 1, }, # use instrument_data.id instead
        model_groups     => { via => 'model', is_many => 1, },      # use model.groups instead
        work_order_names => { via => 'work_orders', to => 'name' }, # use work_orders.name instead
        work_order_numbers => { via => 'work_orders', to => 'id' }, # use work_orders.number instead
        
        # this is just a duplicate of "status" with a different name
        master_event_status     => { via => 'the_master_event', to => 'event_status' },

        # we now use model/build inputs instead of links
        # when these can be removed do
        # see "downstream_builds" 
        from_build_links => { is => 'Genome::Model::Build::Link', reverse_as => 'to_build',
                              doc => 'bridge table entries where this is the \"to\" build(used to retrieve builds this build is \"from\")' },
        from_builds      => { is => 'Genome::Model::Build', via => 'from_build_links', to => 'from_build',
                              doc => 'Genome builds that contribute \"to\" this build' },
        to_build_links   => { is => 'Genome::Model::Build::Link', reverse_as => 'from_build',
                              doc => 'bridge entries where this is the \"from\" build(used to retrieve builds builds this build is \"to\")' },
        to_builds        => { is => 'Genome::Model::Build', via => 'to_build_links', to => 'to_build',
                              doc => 'Genome builds this build contributes \"to\"' },
    ],
    schema_name => 'GMSchema',
    data_source => 'Genome::DataSource::GMSchema',
    id_generator => '-uuid',
};

sub __display_name__ {
    my $self = shift;
    return $self->id . ' of ' . $self->model->name;
}

sub model_class {
    my $self = shift;
    my $model_class = $self->class;
    $model_class =~ s/::Model::Build/::Model/;
    return $model_class;
}

sub data_set_path {
    my ($self, $dataset, $version, $file_format) = @_;
    my $path;
    
    if ($version and $file_format) {
        $version =~ s/^v//;
        $path = $self->data_directory."/$dataset.v$version.$file_format";
    }
    elsif ($file_format) {
        # example $b->data_set_path('alignments/tumor/','','flagstat')
        my @paths = glob($self->data_directory."/$dataset".'*.'.$file_format);
        return unless @paths == 1;
        $path = $paths[0];
    }
    return $path if -e $path;
    return;
}

# TODO Remove this
sub _resolve_subclass_name_by_sequencing_platform { # only temporary, subclass will soon be stored
    my $class = shift;

    Carp::confess("this is used by sub-classes which further subclassify by sequencing platform!")
        if $class eq __PACKAGE__;

    my $sequencing_platform;
    if (ref($_[0]) and $_[0]->isa('Genome::Model::Build')) {
        $sequencing_platform = $_[0]->model->sequencing_platform;
    }
    else {
        my %params;
        if (ref($_[0]) and $_[0]->isa("UR::BoolExpr")) {
            %params = $_[0]->params_list;
        }
        else {
            %params = @_;
        }
        my $model_id = $params{model_id};
        $class->_validate_model_id($params{model_id})
            or return;
        my $model = Genome::Model->get($params{model_id});
        unless ( $model ) {
            Carp::confess("Can't get model for id: .".$params{model_id});
        }
        $sequencing_platform = $model->sequencing_platform;
    }

    return unless $sequencing_platform;

    return $class. '::'.Genome::Utility::Text::string_to_camel_case($sequencing_platform);
}

# auto generate sub-classes for any valid model sub-class
sub __extend_namespace__ {
    # auto generate sub-classes for any valid processing profile
    my ($self,$ext) = @_;

    my $meta = $self->SUPER::__extend_namespace__($ext);
    return $meta if $meta;

    my $model_subclass_name = 'Genome::Model::' . $ext;
    my $model_subclass_meta = eval { $model_subclass_name->__meta__ };
    if ($model_subclass_meta and $model_subclass_name->isa('Genome::Model')) {
        my $build_subclass_name = 'Genome::Model::Build::' . $ext;
        # The actual inputs and metrics are added during subclass definition preprocessing.
        # Then the whole set is expanded, allowing the developer to write less of 
        # # the build class.
        # See Genome/Model/Build.pm _preprocess_subclass_description.
        my $build_subclass_meta = UR::Object::Type->define(
            class_name => $build_subclass_name,
            is => 'Genome::Model::Build',
        );
        die "Error defining $build_subclass_name for $model_subclass_name!" unless $model_subclass_meta;
        return $build_subclass_meta;
    }
    return;
}

sub create {
    my $class = shift;
    if ($class eq __PACKAGE__ or $class->__meta__->is_abstract) {
        # Let the base class re-call the constructor from the correct sub-class
        return $class->SUPER::create(@_);
    }

    my $self = $class->SUPER::create(@_);
    return unless $self;

    eval {
        # Give the model a chance to update itself prior to copying inputs from it
        # TODO: the model can have an observer on it's private subclass of build for
        # the one general case which uses this (Convergence).
        # The other case which uses this is ReferenceAlignment but just for TGI-specific
        # policies, which should really be observers in ::Site::TGI.
        if ($self->model->can("check_for_updates")) {
            unless ($self->model->check_for_updates) {
                Carp::confess "Could not update model!";
            }
        }

        # Now copy (updated) inputs to build
        unless ($self->_copy_model_inputs) {
            Carp::confess "Could not copy model inputs from model " . $self->model->__display_name__ . " to new build!";
        }

        # Allow model to initialize build
        unless ($self->model->_initialize_build($self)) {
            Carp::confess "Model " . $self->model->__display_name__ .
                " could not initialize new build";
        }

        # Create master event, which stores status/user/date created, etc
        unless ($self->_create_master_event) {
            Carp::confess "Could not create master event for new build of model " . $self->model->__display_name__;
        }

        $self->add_note(
            header_text => 'Build Created',
        );
    };

    if ($@) {
        $self->error_message("Could not create new build of model " . $self->__display_name__ . ", reason: $@");
        $self->delete;
        return;
    }

    return $self;
}

sub _create_master_event {
    my $self = shift;
    my $event = Genome::Model::Event->create(
        event_type => 'genome model build',
        event_status => 'New',
        model_id => $self->model->id,
        build_id => $self->id,
    );
    return $event;
}

sub _copy_model_inputs {
    my $self = shift;

    # Failing to copy an input SHOULD NOT be fatal. If the input is required for the build
    # to run, it'll be caught when the build is verified as part of the start method, which
    # will leave the build in an "unstartable" state that can be reviewed later.
    for my $input ($self->model->inputs) {
        eval {
            my %params = map { $_ => $input->$_ } (qw/ name value_class_name value_id filter_desc /);

            # Resolve inputs pointing to a model to a build.
            if($params{value_class_name}->isa('Genome::Model')) {
                my $input_name = $input->name;
                if ($input_name =~ /_model(s)?$/) {
                    $input_name =~ s/_model(?=($|s$))/_build/;
                    $params{name} = $input_name;
                }

                my @existing_inputs = $self->inputs(name => $input_name);
                if (@existing_inputs) {
                    foreach my $existing_input (@existing_inputs) {
                        my $existing_input_value = $existing_input->value;
                        if ($existing_input_value
                            and $existing_input_value->isa('Genome::Model::Build')
                            and $existing_input_value->model_id eq $input->value->id) {
                            die "Input with name $input_name already exists for build!";
                        }
                    }
                }

                my $input_model = $input->value;
                my $input_build = $self->select_build_from_input_model($input_model);
                unless($input_build) {
                    die "Could not resolve a build of model " . $input_model->__display_name__;
                }

                $params{value_class_name} = $input_build->class;
                $params{value_id} = $input_build->id;
            }

            unless ($self->add_input(%params)) {
                die "Could not copy model input " . $params{name} . " with ID " . $params{value_id} .
                    " and class " . $params{value_class_name} . " to new build";
            }
        };
        if ($@) {
            my $error_string = $@;
            $self->warning_message("Could not copy model input " . $input->__display_name__ .
                " to build " . $self->__display_name__ . " of model " . $self->model->__display_name__ .
                " because $error_string");
            next;
        }
    }

    return 1;

}

sub select_build_from_input_model {
    my ($self, $model) = @_;
    return $model->last_complete_build;
}

sub instrument_data_count {
    my $self = shift;
    my @instrument_data = $self->instrument_data;
    if (@instrument_data) {
        return scalar(@instrument_data);
    }
    return 0;
}

# why is this not a defined relationship above? -ss
sub events {
    my $self = shift;
    my @events = Genome::Model::Event->get(
        model_id => $self->model_id,
        build_id => $self->build_id,
        @_,
    );
    return @events;
}

# why is this not a defined relationship above? -ss
sub build_events {
    my $self = shift;
    my @build_events = Genome::Model::Event::Build->get(
        model_id => $self->model_id,
        build_id => $self->build_id,
        @_
    );
    return @build_events;
}

sub build_event {
    my $self = shift;
    my @build_events = $self->build_events;
    if (scalar(@build_events) > 1) {
        my $error_message = 'Found '. scalar(@build_events) .' build events for model id '.
        $self->model_id .' and build id '. $self->build_id ."\n";
        for (@build_events) {
            $error_message .= "\t". $_->desc .' '. $_->event_status ."\n";
        }
        die($error_message);
    }
    return $build_events[0];
}

sub workflow_name {
    my $self = shift;
    return $self->build_id . ' all stages';
}

sub workflow_instances {
    my $self = shift;
    my @instances = Workflow::Operation::Instance->get(
        name => $self->workflow_name,
    );
    return @instances;
}

sub newest_workflow_instance {
    my $self = shift;
    my @sorted = sort {
        $b->id <=> $a->id
    } $self->workflow_instances;
    if (@sorted) {
        return $sorted[0];
    } else {
        return;
    }
}

sub cpu_slot_hours {
    my $self = shift;
    my $breakdown = $self->_cpu_slot_usage_breakdown(@_);
    my $total = 0;
    for my $step (sort keys %$breakdown) {
        my $sum     = $breakdown->{$step}{sum};
        my $count   = $breakdown->{$step}{count};
        print join("\t","BREAKDOWN:",$self->id,$self->model->name,$step,$sum,$count);
        $total += $sum;
    }

    return $total/60;
}

sub _cpu_slot_usage_breakdown {
    my $self = shift;
    my @params = @_;

    my %steps;

    my $workflow_instance = $self->newest_workflow_instance;
    my $bx;
    if (@params) {
        if ($params[0] =~ /event_type/) {
            $params[0] =~ s/event_type/name/;
        }
        $bx = Workflow::Operation::Instance->define_boolexpr(@params);
    }


    my @all_children;
    my @queue = $workflow_instance;

    while ( my $next = shift @queue ) {
        my @children = $next->related_instances;
        push @all_children, @children;
        push @queue,        @children;
    }

    for my $op_inst (@all_children) {
        if ($bx and not $bx->evaluate($op_inst)) {
            warn "skipping $op_inst $op_inst->{name}...\n";
            next;
        }
        my $op = $op_inst->operation;
        my $op_name = $op->name;
        my $op_type = $op->operation_type;

        if ($op_type->class eq 'Workflow::OperationType::Model') {
            # don't double count
            next;
        }

        if ($op_type->can('lsf_queue') and defined($op_type->lsf_queue)
            and $op_type->lsf_queue eq $ENV{GENOME_LSF_QUEUE_BUILD_WORKFLOW}
        ) {
            # skip jobs which run in workflow because they internally run another workflow
            next;
        }

        my $cpus = 1;
        my $rusage = ($op_type->can('lsf_resource') ? $op_type->lsf_resource : undef);
        if (defined($rusage)) {
            if ($rusage =~ /(?<!\S)-n\s+(\S+)/) {
                $cpus = $1;
            }
        }

        my $eclass;
        if ($op_type->can('command_class_name')) {
            $eclass = $op_type->command_class_name;
        }

        my $d1    = Date::Manip::ParseDate( $op_inst->end_time );
        my $d2    = Date::Manip::ParseDate( $op_inst->start_time );
        my $delta = Date::Manip::DateCalc( $d2, $d1 );
        my $value = Date::Manip::Delta_Format( $delta, 1, "%mt" );

        if ($value eq '') {
            $value = 0; # for crashed/incomplete steps
        }

        if ($cpus eq '') {
            die "null cpus??? resource was $rusage\n";
        }

        my $key         = $op_inst->name;
        $key =~ s/ \d+$//;
        $steps{$key}{sum} ||= 0;
        $steps{$key}{count} ||= 0;
        $steps{$key}{sum} += $value * $cpus;
        $steps{$key}{count} += $cpus;
    }

    return \%steps;
}

sub calculate_estimated_kb_usage {
    my $self = shift;

    # Default of 500 MiB in case a subclass fails to
    # override this method.  At least this way there
    # will be an allocation, which will likely be
    # wildly inaccurate, but if the build fails to fail,
    # when it finishes, it will reallocate down to the
    # actual size.  Whereas the previous behaviour
    # (return undef) caused *no* allocation to be made.
    # Which it has been decided is a bigger problem.
    return 512_000;
}

sub get_or_create_data_directory {
    #     If the data directory is not set, resolving it requires making an
    # allocation.  A build is unlikely to make a new allocation at any other
    # time, so a separate build instance method for allocating is not provided.
    my $self = shift;
    return $self->data_directory if $self->data_directory;

    my $allocation_path = 'model_data/' . $self->model->id . '/build'. $self->build_id;
    my $kb_requested = $self->calculate_estimated_kb_usage;
    unless ($kb_requested) {
        $self->error_message("Could not estimate kb usage for allocation!");
        return;
    }

    my $disk_group_name = $self->model->_resolve_disk_group_name_for_build($self);
    unless ($disk_group_name) {
        $self->error_message('Failed to resolve a disk group for a new build!');
        return;
    }

    my $disk_allocation = Genome::Disk::Allocation->create(
        disk_group_name => $disk_group_name,
        allocation_path => $allocation_path,
        kilobytes_requested => $kb_requested,
        owner_class_name => $self->class,
        owner_id => $self->id,
    );
    unless ($disk_allocation) {
        $self->error_message("Could not create allocation for build " . $self->__display_name__);
        return;
    }

    $self->data_directory($disk_allocation->absolute_path);
    return $self->data_directory;
}

sub reallocate {
    my $self = shift;

    my $status = $self->status;
    my $disk_allocation = $self->disk_allocation;

    if ($disk_allocation) {
        my $reallocated = eval { $disk_allocation->reallocate };
        $self->warning_message("Failed to reallocate disk space!") unless $reallocated;
    }
    else {
        $self->warning_message("Reallocate called for build (" . $self->__display_name__ . ") but it does not have a disk allocation.");
    }

    # Always returns 1 due to legacy behavior.
    return 1;
}

sub all_allocations {
    my $self = shift;

    my @allocations;
    push @allocations, $self->disk_allocation if $self->disk_allocation;

    push @allocations, $self->user_allocations;
    push @allocations, $self->symlinked_allocations;

    push @allocations, $self->input_allocations;
    push @allocations, $self->event_allocations;

    my %allocations = map { $_->id => $_ } @allocations;
    return values %allocations;
}

sub input_allocations {
    my $self = shift;

    my @allocations;
    for my $input ($self->inputs) {
        my $value = $input->value;
        if ($value and $value->isa('Genome::Model::Build')) {
            push @allocations, $input->value->all_allocations();
        }
        else {
            push @allocations, Genome::Disk::Allocation->get(
                owner_id => $input->value_id,
                owner_class_name => $input->value_class_name,
            );
        }
    }
    return @allocations;
}

sub user_allocations {
    my $self = shift;

    my @results = $self->all_results;
    my @allocations = map { $_->disk_allocations } @results;
    return @allocations;
}


sub event_allocations {
    my $self = shift;

    my @allocations;
    my @events = $self->events;
    for my $event (@events) {
        my @event_allocations = Genome::Disk::Allocation->get(
            owner_id => $event->id,
            owner_class_name => $event->class,
        );
        push @allocations, @event_allocations if @event_allocations;
    }
    return @allocations;
}


sub all_results {
    my $self = shift;

    my @users = $self->result_users;
    my %seen;
    my @results;
    while(@users) {
        my $u = shift @users;
        my $result = $u->software_result;
        next if $seen{$result->id}; #already processed

        push @results, $result;
        push @users, Genome::SoftwareResult::User->get(user_class_name => $result->class, user_id => $result->id);
        $seen{$result->id}++;
    }

    return @results;
}

sub symlinked_allocations {
    my $self = shift;

    my $data_directory = $self->data_directory;
    return if not $data_directory or not -d $data_directory;

    # In some circumstances, this can get called (though all_allocations())
    # many times.  NFS slowness would compound the problem and make this
    # method way too slow.  We'll get the result the first time and cache
    # the answer
    unless ($self->{__symlinked_allocations}) {

        my %symlinks;
        File::Find::find(
            {
                wanted => sub{
                    return if not -l $File::Find::name;
                    return if -d $File::Find::name;
                    $symlinks{$File::Find::name} = readlink($File::Find::name);
                },
                follow_fast => 1,
                follow_skip => 2,
            },
            $data_directory,
        );

        my %allocations;
        for my $symlink ( keys %symlinks ) {
            my $target = $symlinks{$symlink};
            my @tokens = split(m#/+#, $target);
            my $allocation_path = join('/', @tokens[4..$#tokens]);
            my @allocations = Genome::Disk::Allocation->get(allocation_path => $allocation_path);
            next if not @allocations or @allocations > 1;
            my $allocation = $allocations[0];
            next if $allocations{$allocation->id};
            UR::Context->reload($allocation); # this may have been loaded and unarchived earlier
            $allocations{$allocation->id} = $allocation;
            $allocation->{_symlink_name} = $symlink;
            $allocation->{_target_name} = $target;
            $allocation->{_target_exists} = ( -e $target ? 1 : 0 );
        }

        my @allocs = values %allocations;
        $self->{__symlinked_allocations} = \@allocs;
        $self->{__symlinked_allocations_time} = time();
    }

    return @{ $self->{__symlinked_allocations}};
}

sub input_builds {
    my $self = shift;

    my @builds;
    for my $input ($self->inputs) {
        my $value = $input->value;
        if ($value and $value->isa('Genome::Model::Build')) {
            push @builds, $value;
            push @builds, $value->input_builds;
        }
    }
    return @builds;
}

sub relink_symlinked_allocations {
    my $self = shift;
    $self->debug_message('Relink symlinked allocations...');

    my @symlinked_allocations = $self->symlinked_allocations;
    $self->debug_message('Found '.@symlinked_allocations.' symlinked allocations');
    return 1 if not @symlinked_allocations;

    for my $symlinked_allocation ( @symlinked_allocations ) {
        $self->debug_message('Allocation: '.$symlinked_allocation->id);
        my $symlink_name = $symlinked_allocation->{_symlink_name};
        $self->debug_message("Symlink name: $symlink_name");
        my $target_name = $symlinked_allocation->{_target_name};
        $self->debug_message("Target name: $target_name");
        if ( $symlinked_allocation->{_target_exists} ) {
            $self->debug_message('Target exists! Skipping...');
            next;
        }
        if ( $symlinked_allocation->absolute_path eq $symlinked_allocation->{_target_name} ) {
            $self->debug_message('Target name matches absolute path! Skipping...');
            next;
        }
        my $new_target = $symlinked_allocation->absolute_path;
        if ( $new_target =~ m#^/gscarchive# ) {# preserve link if still archived
            $self->error_message('Allocation failed to reload or is still archived. Please correct. Skipping...');
            next;
        }
        $self->debug_message("Remove broken link: $symlink_name TO $target_name");
        unlink($symlink_name);
        $self->debug_message("Add link: $symlink_name TO $new_target");
        symlink($new_target, $symlink_name);
        $self->error_message('Failed to relink allocation!') if not -l $symlink_name;
    }

    $self->debug_message('Relink symlinked allocations...Done');
    return 1;
}

sub log_directory {
    my $self = shift;
    return unless $self->data_directory;
    return  $self->data_directory . '/logs/';
}

sub reports_directory {
    my $self = shift;
    return unless $self->data_directory;
    return $self->data_directory . '/reports/';
}

sub resolve_reports_directory { return reports_directory(@_); } #????

sub add_report {
    my ($self, $report) = @_;

    my $directory = $self->resolve_reports_directory;
    die "Could not resolve reports directory" unless $directory;
    if (-d $directory) {
        my $subdir = $directory . '/' . $report->name_to_subdirectory($report->name);
        if (-e $subdir) {
            $self->debug_message("Sub-directory $subdir exists!   Moving it out of the way...");
            my $n = 1;
            my $max = 20;
            while ($n < $max and -e $subdir . '.' . $n) {
                $n++;
            }
            if ($n == $max) {
                die "Too many re-runs of this report!  Contact Informatics..."
            }
            rename $subdir, "$subdir.$n";
            if (-e $subdir) {
                die "failed to move old report dir $subdir to $subdir.$n!: $!";
            }
        }
    }
    else {
        $self->status_message("creating directory $directory...");
        unless (Genome::Sys->create_directory($directory)) {
            die "failed to make directory $directory!: $!";
        }
    }

    if ($report->save($directory)) {
        $self->debug_message("Saved report to override directory: $directory");
        return 1;
    }
    else {
        $self->error_message("Error saving report!: " . $report->error_message());
        return;
    }
}

sub archivable {
    my $self = shift;
    my $allocation = $self->disk_allocation;
    unless ($allocation) {
        $self->warning_message("Could not get allocation for build " . $self->__display_name__);
        return 0;
    }
    return $allocation->archivable();
}

sub is_archived {
    my $self = shift;
    my $is_archived = 0;
    my @allocations = $self->all_allocations;
    for my $allocation (@allocations) {
        if ($allocation->is_archived()) {
            $is_archived = 1;
            last;
        }
    }
    return $is_archived;
}

sub start {
    my $self = shift;
    my %params = @_;

    # Regardless of how this goes, build requested should be unset. And we also want to know what software was used.
    $self->model->build_requested(0);
    $self->software_revision(Genome::Sys->snapshot_revision) unless $self->software_revision;

    eval {
        # Validate build for start and collect tags that represent problems.
        # Croak is used here instead of confess to limit error message length. The entire message must fit into the
        # body text of a note, and can cause commit problems if the length exceeds what the db column can accept.
        # TODO Delegate to some other method to create the error message
        my @tags = $self->validate_for_start;
        if (@tags) {
            my @msgs;
            for my $tag (@tags) {
                push @msgs, $tag->__display_name__;
            }
            Carp::croak "Build " . $self->__display_name__ . " could not be validated for start!\n" . join("\n", @msgs);
        }

        # Either returns the already-set data directory or creates an allocation for the data directory
        unless ($self->get_or_create_data_directory) {
            Carp::croak "Build " . $self->__display_name__ . " failed to resolve a data directory!";
        }

        # Give builds an opportunity to do some initialization after the data directory has been resolved
        unless ($self->post_allocation_initialization) {
            Carp::croak "Build " . $self->__display_name__ . " failed to initialize after resolving data directory!";
        }

        # Record that the event has been scheduled.
        $self->the_master_event->schedule;

        # Creates a workflow for the build
        # TODO Initialize workflow shouldn't take arguments
        unless ($self->_initialize_workflow($params{job_dispatch} || $ENV{GENOME_LSF_QUEUE_BUILD_WORKER_ALT})) {
            Carp::croak "Build " . $self->__display_name__ . " could not initialize workflow!";
        }

        # Launches the workflow (in a pend state, it's resumed by a commit hook)
        unless ($self->_launch(%params)) {
            Carp::croak "Build " . $self->__display_name__ . " could not be launched!";
        }

        $self->add_note(
            header_text => 'Build Started',
        );
    };

    if ($@) {
        my $error = $@;
        $self->add_note(
            header_text => 'Unstartable',
            body_text => "Could not start build, reason: $error",
        );
        $self->the_master_event->event_status('Unstartable');
        $self->error_message("Could not start build " . $self->__display_name__ . ", reason: $error");
        return;
    }

    return 1;
}

sub post_allocation_initialization {
    #     Override in subclasses to allow build to do initialization just after
    # 'data_directory' has been allocated.
    return 1;
}

sub validate_for_start_methods {
    # Be very wary of removing any of these as many subclasses use SUPER::validate_for_start_methods.
    # Each method should return tags.
    my @methods = (
        # 'validate_inputs_have_values' should be checked first.
        'validate_inputs_have_values',
        'inputs_have_compatible_reference',
        'validate_instrument_data',
        # instrument_data_assigned, # Several build subclasses use this.
    );
    return @methods;
}

sub validate_for_start {
    #     Run  all the 'validate_for_start_methods' and collect tags that
    # correspond to errors.
    my $self = shift;

    my @tags;
    my @methods = $self->validate_for_start_methods;

    for my $method (@methods) {
        unless ($self->can($method)) {
            die $self->warning_message("Validation method $method not found!");
        }
        my @returned_tags = grep { defined $_ } $self->$method(); # Prevents undef from being pushed to tags list
        push @tags, @returned_tags if @returned_tags;
    }

    return @tags;
}

sub instrument_data_assigned {
    #     Since this could be used by several build subclasses it is here
    # but it is not a default 'validate_for_start_method' for all builds.
    my $self = shift;
    my @tags;
    my @instrument_data = $self->instrument_data;
    unless (@instrument_data) {
        push @tags, UR::Object::Tag->create(
            type => 'error',
            properties => ['instrument_data'],
            desc => 'no instrument data assigned to build',
        );
    }
    return @tags;
}

sub validate_instrument_data{
    my $self = shift;
    my @tags;
    my @instrument_data = $self->instrument_data;
    for my $instrument_data (@instrument_data){
        if (not defined $instrument_data->read_count){
            push @tags, UR::Object::Tag->create(
                type => 'error',
                properties => ['instrument_data'],
                desc => 'read count for instrument data (' . $instrument_data->id . ') has not been calculated, use `genome instrument-data calculate-read-count` to correct this.',
            );
        } elsif (not $instrument_data->read_count){
            push @tags, UR::Object::Tag->create(
                type => 'error',
                properties => ['instrument_data'],
                desc => 'no reads for instrument data (' . $instrument_data->id . ') assigned to build',
            );
        }
    }
    return @tags;
}

sub validate_inputs_have_values {
    my $self = shift;
    my @inputs = grep { $self->can($_->name) } $self->inputs;

    # find inputs with undefined values
    my @inputs_without_values = grep { not defined $_->value } @inputs;
    my %input_names_to_ids;
    for my $input (@inputs_without_values){
        my $name = $input->name;
        if ($self->can("name")) {
            my $pmeta = $self->__meta__->property($name);
            if ($pmeta) {
                if ($pmeta->is_optional) {
                    next;
                }
            }
        }
        $input->value;
        $input_names_to_ids{$input->name} .= $input->value_class_name . ":" . $input->value_id . ',';
    }

    # make tags for inputs with undefined values
    my @tags;
    for my $input_name (keys %input_names_to_ids) {
        push @tags, UR::Object::Tag->create(
            type => 'error',
            properties => [$input_name],
            desc => "Value no longer exists for value id: " . $input_names_to_ids{$input_name},
        );
    }

    return @tags;
}

sub inputs_have_compatible_reference {
    my $self = shift;

    my @reference_sequence_methods = ('reference_sequence', 'reference', 'reference_sequence_build');

    # Determine the reference sequence for the build.
    my ($build_reference_method) = grep { $self->can($_) } @reference_sequence_methods;
    return unless $build_reference_method;
    my $build_reference_sequence = $self->$build_reference_method;

    # Determine the reference sequence for the inputs (and ensure compatiblity).
    my @inputs = $self->inputs;
    my @incompatible_properties;
    for my $input (@inputs) {
        my $object = $input->value;
        next unless $object; #this is reported in validate_inputs_have_values
        my ($input_reference_method) = grep { $object->can($_) } @reference_sequence_methods;
        next unless $input_reference_method;
        my $object_reference_sequence = $object->$input_reference_method;
        next unless $object_reference_sequence;

        unless($object_reference_sequence->is_compatible_with($build_reference_sequence)
                or $self->reference_being_replaced_for_input($input)) {
            push @incompatible_properties, $input->name;
        }
    }

    # Create tags for all the inputs with incompatible reference sequences.
    my $tag;
    if (@incompatible_properties) {
        $tag = UR::Object::Tag->create(
            type => 'error',
            properties => \@incompatible_properties,
            desc => "Not compatible with build's reference sequence '" . $build_reference_sequence->__display_name__ . "'.",
        );
    }

    return $tag;
}

sub reference_being_replaced_for_input {
    my $self = shift;
    my $input = shift;

    #for overriding in subclasses--by default none are replaced
    #(example of when this would be true: an imported BAM being realigned)
    return;
}


sub stop {
    my $self = shift;

    $self->status_message('Attempting to stop build: '.$self->id);

    my $user = getpwuid($<);
    if ($user ne 'apipe-builder' && $user ne $self->run_by) {
        $self->error_message("Can't stop a build originally started by: " . $self->run_by);
        return 0;
    }

    my $job = $self->_get_running_master_lsf_job;
    if ( defined $job ) {
        $self->status_message('Killing job: '.$job->{Job});
        $self->_kill_job($job);
        $self = Genome::Model::Build->load($self->id);
    }

    $self->add_note(
        header_text => 'Build Stopped',
    );

    my $self_event = $self->build_event;
    my $error = Genome::Model::Build::Error->create(
        build_event_id => $self_event->id,
        stage_event_id => $self_event->id,
        stage => 'all stages',
        step_event_id => $self_event->id,
        step => 'main',
        error => 'Killed by user',
    );

    $self->status_message('Failing build: '.$self->id);
    unless ($self->fail($error)) {
        $self->error_message('Failed to fail build');
        return;
    }

    return 1
}

sub _kill_job {
    my ($self, $job) = @_;

    Genome::Sys->shellcmd(
        cmd => 'bkill '.$job->{Job},
    );

    my $i = 0;
    do {
        $self->status_message("Waiting for job to stop") if ($i % 10 == 0);
        $i++;
        sleep 1;
        $job = $self->_get_job( $job->{Job} );

        if ($i > 60) {
            $self->error_message("Build master job did not die after 60 seconds.");
            return 0;
        }
    } while ($job && ($job->{Status} ne 'EXIT' && $job->{Status} ne 'DONE'));

    return 1;
}

sub _get_running_master_lsf_job {
    my $self = shift;

    my $job_id = $self->the_master_event->lsf_job_id;
    return if not defined $job_id;

    my $job = $self->_get_job($job_id);
    return if not defined $job;

    if ( $job->{Status} eq 'EXIT' or $job->{Status} eq 'DONE' ) {
        return;
    }

    return $job;
}

sub _get_job {
    use Genome::Model::Command::Services::Build::Scan;
    my $self = shift;
    my $job_id = shift;

    my @jobs = ();
    my $iter = Job::Iterator->new($job_id);
    while (my $job = $iter->next) {
        push @jobs, $job;
    }

    if (@jobs > 1) {
        $self->error_message("More than 1 job found for this build? Alert apipe");
        return 0;
    }

    return shift @jobs;
}

sub _launch {
    my $self = shift;
    my %params = @_;

    local $ENV{UR_DUMP_DEBUG_MESSAGES} = 1;
    local $ENV{UR_COMMAND_DUMP_DEBUG_MESSAGES} = 1;
    local $ENV{UR_DUMP_STATUS_MESSAGES} = 1;
    local $ENV{UR_COMMAND_DUMP_STATUS_MESSAGES} = 1;
    local $ENV{GENOME_BUILD_ID} = $self->id;

    # right now it is "inline" or the name of an LSF queue.
    # ultimately, it will be the specification for parallelization
    # including whether the server is inline, forked, or bsubbed, and the
    # jobs are inline, forked or bsubbed from the server
    my $model = $self->model;

    my $server_dispatch = _server_dispatch($model, \%params);
    my $job_dispatch = _job_dispatch($model, \%params);
    my $job_group_spec = _job_group_spec(\%params);

    # all params should have been deleted (as they were handled)
    die "Bad params!  Expected server_dispatch and job_dispatch!" . Data::Dumper::Dumper(\%params) if %params;

    my $build_event = $self->the_master_event;

    if ($server_dispatch eq 'inline') {
        my %args = (
            model_id => $self->model_id,
            build_id => $self->id,
        );
        if ($job_dispatch eq 'inline') {
            $args{inline} = 1;
        }

        my $rv = Genome::Model::Command::Services::Build::Run->execute(%args);
        return $rv;
    }
    else {
        my $add_args = ($job_dispatch eq 'inline') ? ' --inline' : '';

        # bsub into the queue specified by the dispatch spec
        my $lsf_project = "build" . $self->id;
        $ENV{'WF_LSF_PROJECT'} = $lsf_project;
        my $user = Genome::Sys->username;
        my $lsf_command  = join(' ',
            'bsub -N -H',
            '-P', $lsf_project,
            '-q', $server_dispatch,
            ($ENV{WF_EXCLUDE_JOB_GROUP} ? '' : $job_group_spec),
            '-u', Genome::Utility::Email::construct_address(),
            '-o', $build_event->output_log_file,
            '-e', $build_event->error_log_file,
            'annotate-log',
            ($Command::entry_point_bin || 'genome'),
            'model services build run',
            $add_args,
            '--model-id', $model->id,
            '--build-id', $self->id,
        );
        my $job_id = $self->_execute_bsub_command($lsf_command);
        return unless $job_id;

        $build_event->lsf_job_id($job_id);

        return 1;
    }
}

sub _job_dispatch {
    my $model = shift;
    my $params = shift;
    my $job_dispatch;
    if (exists($params->{job_dispatch})) {
        $job_dispatch = delete $params->{job_dispatch};
    } elsif ($model->processing_profile->can('job_dispatch') && defined $model->processing_profile->job_dispatch) {
        $job_dispatch = $model->processing_profile->job_dispatch;
    } else {
        $job_dispatch = $ENV{GENOME_LSF_QUEUE_BUILD_WORKER_ALT};
    }
    return $job_dispatch;
}

sub _server_dispatch {
    my $model = shift;
    my $params = shift;
    my $server_dispatch;
    if (exists($params->{server_dispatch})) {
        $server_dispatch = delete $params->{server_dispatch};
    } elsif ($model->processing_profile->can('server_dispatch') && defined $model->processing_profile->server_dispatch) {
        $server_dispatch = $model->processing_profile->server_dispatch;
    } elsif ($model->can('server_dispatch') && defined $model->server_dispatch) {
        $server_dispatch = $model->server_dispatch;
    } else {
        $server_dispatch = $ENV{GENOME_LSF_QUEUE_BUILD_WORKFLOW};
    }
    return $server_dispatch;
}

sub _default_job_group {
    my $user = getpwuid($<);
    return '/apipe-build/' . $user;
}

sub _job_group_spec {
    my $params = shift;
    my $job_group = _default_job_group();
    if (exists $params->{job_group}) {
        $job_group = delete $params->{job_group};
    }
    return ($job_group ? " -g $job_group" : '');
}

sub _initialize_workflow {
    #     Create the data and log directories and resolve the workflow for this build.
    my $self = shift;
    my $optional_lsf_queue = shift || $ENV{GENOME_LSF_QUEUE_BUILD_WORKER_ALT};

    Genome::Sys->create_directory( $self->data_directory )
        or return;

    Genome::Sys->create_directory( $self->log_directory )
        or return;

    my $model = $self->model;
    my $processing_profile = $self->processing_profile;
    my $workflow = $model->_resolve_workflow_for_build($self, $optional_lsf_queue);

    ## so developers dont fail before the workflow changes get deployed to /gsc/scripts
    # NOTE: Genome::Config is obsolete, so this code must work when it is not installed as well.
    if ($workflow->can('notify_url') and $ENV{GENOME_SYS_SERVICES_WEB_VIEW_URL}) {
        require UR::Object::View::Default::Xsl;

        my $cachetrigger = $ENV{GENOME_SYS_SERVICES_WEB_VIEW_URL};
        $cachetrigger =~ s/view$/cachetrigger/;

        my $url = $cachetrigger . '/' . UR::Object::View::Default::Xsl::type_to_url(ref($self)) . '/status.html?id=' . $self->id;
        $url .= ' ' . $cachetrigger . '/workflow/operation/instance/statuspopup.html?id=[WORKFLOW_ID]';

        $workflow->notify_url($url);
    }
    $workflow->save_to_xml(OutputFile => $self->data_directory . '/build.xml');

    return $workflow;
}

sub _execute_bsub_command { # here to overload in testing
    my ($self, $cmd) = @_;

    local $ENV{UR_DUMP_DEBUG_MESSAGES} = 1;
    local $ENV{UR_COMMAND_DUMP_DEBUG_MESSAGES} = 1;
    local $ENV{UR_DUMP_STATUS_MESSAGES} = 1;
    local $ENV{UR_COMMAND_DUMP_STATUS_MESSAGES} = 1;

    if ($ENV{UR_DBI_NO_COMMIT}) {
        $self->warning_message("Skipping bsub when NO_COMMIT is turned on (job will fail)\n$cmd");
        return 1;
    }

    my $bsub_output = `$cmd`;

    my $rv = $? >> 8;
    if ( $rv ) {
        $self->error_message("Failed to launch bsub (exit code: $rv) command:\n$bsub_output");
        return;
    }

    if ( $bsub_output =~ m/Job <(\d+)>/ ) {
        my $job_id = $1;

        # create a change record so that if it is "undone" it will kill the job
        my $bsub_undo = sub {
            $self->status_message("Killing LSF job ($job_id) for build " . $self->__display_name__ . ".");
            system("bkill $job_id");
        };
        my $lsf_change = UR::Context::Transaction->log_change($self, 'UR::Value', $job_id, 'external_change', $bsub_undo);
        unless ($lsf_change) {
            die $self->error_message("Failed to record LSF job submission ($job_id).");
        }

        # create a commit observer to resume the job when build is committed to database
        my $process = UR::Context->process;
        my $commit_observer = $process->add_observer(
            aspect => 'commit',
            callback => sub {
                my $bresume_output = `bresume $job_id`; chomp $bresume_output;
                $self->status_message($bresume_output) unless ( $bresume_output =~ /^Job <$job_id> is being resumed$/ );
            },
        );
        unless ($commit_observer) {
            $self->error_message("Failed to add commit observer to resume LSF job ($job_id).");
        }

        return "$job_id";
    }
    else {
        $self->error_message("Launched busb command, but unable to parse bsub output: $bsub_output");
        return;
    }
}

sub initialize {
    my $self = shift;

    $self->_verify_build_is_not_abandoned_and_set_status_to('Running')
        or return;

    $self->generate_send_and_save_report('Genome::Model::Report::BuildInitialized')
        or return;

    return 1;
}

sub fail {
    my ($self, @errors) = @_;

    # reload all the events
    my @e = Genome::Model::Event->load(build_id => $self->build_id);

    $self->_verify_build_is_not_abandoned_and_set_status_to('Failed', 1)
        or return;

    if ($self->disk_allocation) {
        $self->reallocate;
    }

    # set event status
    for my $e ($self->the_events(event_status => 'Running')) {
        $e->event_status('Failed');
    }

    $self->generate_send_and_save_report(
        'Genome::Model::Report::BuildFailed', {
            errors => \@errors,
        },
    )
        or return;

    for my $error (@errors) {
        $self->add_note(
            header_text => 'Failed Stage',
            body_text => $error->stage,
        );
        $self->add_note(
            header_text => 'Failed Step',
            body_text => $error->step,
        );
        $self->add_note(
            header_text => 'Failed Error',
            body_text => $error->error,
            auto_truncate_body_text => 1,
        );
    }

    return 1;
}

sub success {
    my $self = shift;

    # reload all the events
    my @e = Genome::Model::Event->load(build_id => $self->build_id);

    # set status
    $self->_verify_build_is_not_abandoned_and_set_status_to('Succeeded', 1)
        or return;

    # set event status
    for my $e ($self->the_events(event_status => ['Running','Scheduled'])) {
        $e->event_status('Abandoned');
    }

    # report - if this fails set status back to Running, then the workflow will fail it
    unless ( $self->generate_send_and_save_report( $self->report_generator_class_for_success ) ) {
        $self->_verify_build_is_not_abandoned_and_set_status_to('Running');
        return;
    }

    # Launch new builds for any convergence models containing this model.
    # To prevent infinite loops, don't do this for convergence builds.
    # FIXME convert this to use the commit callback and model links with a custom notify that doesn't require succeeded builds
    if($self->type_name !~ /convergence/) {
        for my $model_group ($self->model->model_groups) {
            eval {
                $model_group->schedule_convergence_rebuild;
            };
            if($@) {
                $self->error_message('Could not schedule convergence build for model group ' . $model_group->id . '.  Continuing anyway.');
            }
        }
    }

    my $commit_callback;
    $commit_callback = sub {
        $self->the_master_event->cancel_change_subscription('commit', $commit_callback); #only fire once
        $self->debug_message('Firing build success commit callback.');
        my $result = eval {
            Genome::Search->queue_for_update($self->model);
            $self->model->_trigger_downstream_builds($self);
        };
        if($@) {
            $self->error_message('Error executing success callback: ' . $@);
            return;
        }
        unless($result) {
            $self->error_message('Success callback failed.');
            return;
        }

        return UR::Context->commit; #a separate commit is necessary for any changes in the callback
    };

    #The build itself has no __changes__ and UR::Context->commit() will not trigger the subscription if on that object, so
    #use the master build event which has just been updated to 'Succeeded' with the current time.
    $self->the_master_event->create_subscription(
        method => 'commit',
        callback => $commit_callback,
    );

    # reallocate - always returns true (legacy behavior)
    $self->reallocate;

    # TODO Reconsider this method name
    $self->perform_post_success_actions;

    return 1;
}

# TODO Reconsider this name
sub perform_post_success_actions {
    my $self = shift;
    return 1;
}

sub _verify_build_is_not_abandoned_and_set_status_to {
    my ($self, $status, $set_date_completed) = @_;

    my $build_event = $self->build_event;
    # Do we have a master event?
    unless ( $build_event ) {
        $self->error_message(
            'Cannot set build ('.$self->id.") status to '$status' because it does not have a master event."
        );
        return;
    }

    # Is it abandoned?
    if ( $build_event->event_status eq 'Abandoned' ) {
        $self->error_message(
            'Cannot set build ('.$self->id.") status to '$status' because the master event has been abandoned."
        );
        return;
    }

    # Set status and date completed
    $build_event->event_status($status);
    $build_event->date_completed( UR::Context->current->now ) if $set_date_completed;

    return $build_event;
}


sub abandon {
    my $self = shift;
    my $header_text = shift || 'Build Abandoned';
    my $body_text = shift;

    my $status = $self->status;
    if ($status && $status eq 'Abandoned') {
        return 1;
    }

    if ($status && ($status eq 'Running' || $status eq 'Scheduled')) {
        $self->stop;
    }

    # Abandon events
    $self->_abandon_events
        or return;

    # Reallocate - always returns true (legacy behavior)
    if ($self->disk_allocation) {
        $self->reallocate;
    }

    $self->_deactivate_software_results;

    my %add_note_args = (header_text => $header_text);
    $add_note_args{body_text} = $body_text if defined $body_text;
    $self->add_note(%add_note_args);

    Genome::Search->queue_for_update($self->model);

    return 1;
}

sub _deactivate_software_results {
    my $self = shift;
    map{$_->active(0)} Genome::SoftwareResult::User->get(user_class_name => $self->subclass_name, user_id => $self->id);
    return 1;
}

sub _unregister_software_results {
    my $self = shift;
    my @registrations = Genome::SoftwareResult::User->get(user_class_name => $self->subclass_name, user_id => $self->id);
    for my $registration (@registrations){
        unless($registration->delete){
            $self->error_message("Failed to delete registration: " . Data::Dumper::Dumper($registration));
            return;
        }
    }

    return 1;
}

sub _abandon_events { # does not realloc
    my $self = shift;

    my @events = do {
        no warnings;
        sort { $b->date_scheduled cmp $a->date_scheduled } $self->events;
    };

    for my $event ( @events ) {
        unless ( $event->abandon ) {
            $self->error_message(
                sprintf(
                    'Failed to abandon build (%s) because could not abandon event (%s).',
                    $self->id,
                    $event->id,
                )
            );
            return;
        }
    }

    return 1;
}

sub reports {
    my $self = shift;
    my $report_dir = $self->resolve_reports_directory;
    return unless -d $report_dir;
    return Genome::Report->create_reports_from_parent_directory($report_dir);
}

sub get_report {
    my ($self, $report_name) = @_;

    unless ( $report_name ) { # die?
        $self->error_message("No report name given to get report");
        return;
    }

    my $report_dir = $self->reports_directory.'/'.
    Genome::Report->name_to_subdirectory($report_name);
    return unless -d $report_dir;

    return Genome::Report->create_report_from_directory($report_dir);
}

sub available_reports {
    my $self = shift;
    my $report_dir = $self->resolve_reports_directory;
    return unless -d $report_dir;
    return Genome::Report->create_reports_from_parent_directory($report_dir);
}

sub generate_send_and_save_report {
    my ($self, $generator_class, $additional_params) = @_;

    $additional_params ||= {};
    my $generator = $generator_class->create(
        build_id => $self->id,
        %$additional_params,
    );
    unless ( $generator ) {
        $self->error_message(
            sprintf(
                "Can't create report generator (%s) for build (%s)",
                $generator_class,
                $self->id
            )
        );
        return;
    }

    my $report = $generator->generate_report;
    unless ( $report ) {
        $self->error_message(
            sprintf("Can't generate report (%s) for build (%s)", $generator->name, $self->id)
        );
        return;
    }
    $self->add_report($report)
        or return;

    # Do not send report if user is apipe-builder/tester and it is a init, fail or succ report
    my $username = $self->build_event->user_name;
    return 1 if ( $username eq 'apipe-builder' or $username eq 'apipe-tester' )
        and grep { $generator_class eq 'Genome::Model::Report::'.$_ } (qw/ BuildInitialized BuildSucceeded BuildFailed /);
    # or user does not exist
    return 1 if not getpwnam($username);

    my $to = Genome::Sys::User->get(username => $username);
    return 1 if not $to;

    my $email_confirmation = Genome::Report::Email->send_report(
        report => $report,
        to => $to->id,
        from => 'apipe@'.Genome::Config::domain(),
        replyto => 'noreply@'.Genome::Config::domain(),
        # maybe not the best/correct place for this information but....
        xsl_files => [ $generator->get_xsl_file_for_html ],
    );
    unless ( $email_confirmation ) {
        $self->error_message('Couldn\'t email build report ('.lc($report->name).')');
        return;
    }

    return $report;
}

sub report_generator_class_for_success { # in subclass replace w/ summary or the like?
    return 'Genome::Model::Report::BuildSucceeded';
}

#< SUBCLASSING >#
#
# This is called by the infrastructure to appropriately classify abstract processing profiles
# according to their type name because of the "sub_classification_method_name" setting
# in the class definiton...
sub _resolve_subclass_name {
    my $class = shift;

    my $type_name;
	if ( ref($_[0]) and $_[0]->isa(__PACKAGE__) ) {
		$type_name = $_[0]->model->type_name;
	}
    else {
        my ($bx,@extra) = $class->define_boolexpr(@_);
        my %params = ($bx->params_list, @extra);
        my $model_id = $params{model_id};
        my $model = Genome::Model->get($model_id);
        unless ($model) {
            return undef;
        }
        $type_name = $model->type_name;
    }

    unless ( $type_name ) {
        my $rule = $class->define_boolexpr(@_);
        $type_name = $rule->value_for('type_name');
    }

    if (defined $type_name ) {
        my $subclass_name = $class->_resolve_subclass_name_for_type_name($type_name);
        my $sub_classification_method_name = $class->__meta__->sub_classification_method_name;
        if ( $sub_classification_method_name ) {
            if ( $subclass_name->can($sub_classification_method_name)
                 eq $class->can($sub_classification_method_name) ) {
                return $subclass_name;
            } else {
                return $subclass_name->$sub_classification_method_name(@_);
            }
        } else {
            return $subclass_name;
        }
    } else {
        return undef;
    }
}

sub _resolve_subclass_name_for_type_name {
    my ($class,$type_name) = @_;
    my @type_parts = split(' ',$type_name);

    my @sub_parts = map { ucfirst } @type_parts;
    my $subclass = join('',@sub_parts);

    my $class_name = join('::', 'Genome::Model::Build' , $subclass);
    return $class_name;

}

sub _resolve_type_name_for_class {
    my $class = shift;

    my ($subclass) = $class =~ /^Genome::Model::Build::([\w\d]+)$/;
    return unless $subclass;

    return lc join(" ", ($subclass =~ /[a-z\d]+|[A-Z\d](?:[A-Z\d]+|[a-z]*)(?=$|[A-Z\d])/gx));

    my @words = $subclass =~ /[a-z\d]+|[A-Z\d](?:[A-Z\d]+|[a-z]*)(?=$|[A-Z\d])/gx;
    return lc(join(" ", @words));
}

sub get_all_objects {
    my $self = shift;

    my $sorter = sub { # not sure why we sort, but I put it in a anon sub for convenience
        return unless @_;
        #if ( $_[0]->id =~ /^\-/) {
            return sort {$b->id cmp $a->id} @_;
            #}
            #else {
            #return sort {$a->id cmp $b->id} @_;
            #}
    };

    return map { $sorter->( $self->$_ ) } (qw(events inputs metrics from_build_links to_build_links));
}

sub yaml_string {
    my $self = shift;
    my $string = YAML::Dump($self);
    for my $object ($self->get_all_objects) {
        $string .= YAML::Dump($object);
    }
    return $string;
}

sub add_to_build{
    my $self = shift;
    my (%params) = @_;
    my $build = delete $params{to_build};
    my $role = delete $params{role};
    $role||='member';

    $self->error_message("no to_build provided!") and die unless $build;
    my $from_id = $self->id;
    my $to_id = $build->id;
    unless( $to_id and $from_id){
        $self->error_message ( "no value for this build(from_build) id: <$from_id> or to_build id: <$to_id>");
        die;
    }
    my $reverse_bridge = Genome::Model::Build::Link->get(from_build_id => $to_id, to_build_id => $from_id);
    if ($reverse_bridge){
        my $string =  "A build link already exists for these two builds, and in the opposite direction than you specified:\n";
        $string .= "to_build: ".$reverse_bridge->to_build." (this build)\n";
        $string .= "from_build: ".$reverse_bridge->from_build." (the build you are trying to set as a 'to' build for this one)\n";
        $string .= "role: ".$reverse_bridge->role;
        $self->error_message($string);
        die;
    }
    my $bridge = Genome::Model::Build::Link->get(from_build_id => $from_id, to_build_id => $to_id);
    if ($bridge){
        my $string =  "A build link already exists for these two builds:\n";
        $string .= "to_build: ".$bridge->to_build." (the build you are trying to set as a 'to' build for this one)\n";
        $string .= "from_build: ".$bridge->from_build." (this build)\n";
        $string .= "role: ".$bridge->role;
        $self->error_message($string);
        die;
    }
    $bridge = Genome::Model::Build::Link->create(from_build_id => $from_id, to_build_id => $to_id, role => $role);
    return $bridge;
}

sub add_from_build { # rename "add an underlying build" or something...
    my $self = shift;
    my (%params) = @_;
    my $build = delete $params{from_build};
    my $role = delete $params{role};
    $role||='member';

    $self->error_message("no from_build provided!") and die unless $build;
    my $to_id = $self->id;
    my $from_id = $build->id;
    unless( $to_id and $from_id){
        $self->error_message ( "no value for this build(to_build) id: <$to_id> or from_build id: <$from_id>");
        die;
    }
    my $reverse_bridge = Genome::Model::Build::Link->get(from_build_id => $to_id, to_build_id => $from_id);
    if ($reverse_bridge){
        my $string =  "A build link already exists for these two builds, and in the opposite direction than you specified:\n";
        $string .= "to_build: ".$reverse_bridge->to_build." (the build you are trying to set as a 'from' build for this one)\n";
        $string .= "from_build: ".$reverse_bridge->from_build." (this build)\n";
        $string .= "role: ".$reverse_bridge->role;
        $self->error_message($string);
        die;
    }
    my $bridge = Genome::Model::Build::Link->get(from_build_id => $from_id, to_build_id => $to_id);
    if ($bridge){
        my $string =  "A build link already exists for these two builds:\n";
        $string .= "to_build: ".$bridge->to_build." (this build)\n";
        $string .= "from_build: ".$bridge->from_build." (the build you are trying to set as a 'from' build for this one)\n";
        $string .= "role: ".$bridge->role;
        $self->error_message($string);
        die;
    }
    $bridge = Genome::Model::Build::Link->create(from_build_id => $from_id, to_build_id => $to_id, role => $role);
    return $bridge;
}

sub delete {
    my $self = shift;

    # Abandon events
    $self->status_message("Abandoning events associated with build");
    unless ($self->_abandon_events) {
        $self->error_message(
            "Unable to delete build (".$self->id.") because the events could not be abandoned"
        );
        confess $self->error_message;
    }

    # Delete all associated objects
    $self->status_message("Deleting other objects associated with build");
    my @objects = $self->get_all_objects; # TODO this method name should be changed
    for my $object (@objects) {
        $object->delete;
    }

    # Remove the build as a Software Result User
    $self->status_message("Unregistering software results associated with build");
    $self->_unregister_software_results;

    # Deallocate build directory, which will also remove it (unless no commit is on)
    my $disk_allocation = $self->disk_allocation;
    if ($disk_allocation) {
        $self->status_message("Deallocating build directory");
        unless ($disk_allocation->deallocate) {
            $self->warning_message('Failed to deallocate disk space.');
        }
    }

    return $self->SUPER::delete;
}

sub set_metric {
    my $self = shift;
    my $metric_name  = shift;
    my $metric_value = shift;

    my $metric = Genome::Model::Metric->get(build_id=>$self->id, name=>$metric_name);
    my $new_metric;
    if ($metric) {
        #delete an existing one and create the new one
        $metric->delete;
        $new_metric = Genome::Model::Metric->create(build_id=>$self->id, name=>$metric_name, value=>$metric_value);
    } else {
        $new_metric = Genome::Model::Metric->create(build_id=>$self->id, name=>$metric_name, value=>$metric_value);
    }

    return $new_metric->value;
}

sub get_metric {
    my $self = shift;
    my $metric_name = shift;

    my $metric = Genome::Model::Metric->get(build_id=>$self->id, name=>$metric_name);
    if ($metric) {
        return $metric->value;
    }
}

# Returns a list of files contained in the build's data directory
sub files_in_data_directory {
    my $self = shift;
    my @files;
    my $iter = File::Next::files($self->data_directory);
    while(defined (my $file = $iter->())) {
        push @files, $file;
    }
    return \@files;
}

# Given a full path to a file, return a path relative to the build directory
sub full_path_to_relative {
    my ($self, $path) = @_;
    my $rel_path = $path;
    my $dir = $self->data_directory;
    $dir .= '/' unless substr($dir, -1, 1) eq '/';
    $rel_path =~ s/$dir//;
    $rel_path .= '/' if -d $path and substr($rel_path, -1, 1) ne '/';
    return $rel_path;
}

#
# Diff Methods
#
# There has been a slow migration away from individually written build subclasses.
# Ony one of these methods has been updated to defer to the model definition.
# The others will need some sort of updating.
#

# Returns a list of files that should be ignored by the diffing done by compare_output
# Files should be relative to the data directory of the build and can contain regex.
# Override in subclasses!
sub files_ignored_by_diff {
    my $self = shift;
    return $self->model_class->files_ignored_by_build_diff($self);
}

# Returns a list of directories that should be ignored by the diffing done by compare_output
# Directories should be relative to the data directory of the build and can contain regex.
# Override in subclasses!
sub dirs_ignored_by_diff {
    return ();
}

# A list of regexes that, when applied to file paths that are relative to the build's data
# directory, return only one result. This is useful for files that don't have consistent
# names between builds (for example, if they have the build_id embedded in them. Override
# in subclasses!
sub regex_files_for_diff {
    return ();
}

# A list of metrics that the differ should ignore. Some model/build types store information
# as metrics that need to be diffed. Override this in subclasses.
sub metrics_ignored_by_diff {
    return ();
}

# A hash of method suffixes and a file name regex that triggers a custom diff method. This should include those
# files that have timestamps or other changing fields in them that an md5sum can't handle.
# Each suffix should have a method called diff_<SUFFIX> that'll contain the logic.
sub regex_for_custom_diff {
    my $self = shift;

    # Standard custom differs
    my @regex_for_custom_diff = (
        hq     => '\.hq$',
        gz     => '(?<!\.vcf)\.gz$',
        vcf    => '\.vcf$',
        vcf_gz => '\.vcf\.gz$',
    );

    # Addition custom differs
    my $model_class = $self->model_class;
    if ( $model_class->can('addtional_regex_for_custom_diff') ) {
        my @addtional_regex_for_custom_diff = $model_class->addtional_regex_for_custom_diff;
        if ( @addtional_regex_for_custom_diff % 2 != 0 ) {
            Carp::confess('Invalid addtional_regex_for_custom_diff! '.Data::Dumper::Dumper(\@addtional_regex_for_custom_diff));
        }
        push @regex_for_custom_diff, @addtional_regex_for_custom_diff;
    }

    return @regex_for_custom_diff;
}

sub matching_regex_for_custom_diff {
    my $self = shift;
    my $path = shift;

    my %regex_for_custom_diff = $self->regex_for_custom_diff;
    my %matching_regex_for_custom_diff;
    for my $key (keys %regex_for_custom_diff) {
        my $regex = $regex_for_custom_diff{$key};
        $matching_regex_for_custom_diff{$key} = $regex if $path =~ /$regex/;
    }

    return %matching_regex_for_custom_diff;
}

# Gzipped files contain the timestamp and name of the original file, so this prints
# the uncompressed file to STDOUT and pipes it to md5sum.
sub diff_gz {
    my ($self, $first_file, $second_file) = @_;
    my $first_md5  = `gzip -dc $first_file | md5sum`;
    my $second_md5 = `gzip -dc $second_file | md5sum`;
    return 1 if $first_md5 eq $second_md5;
    return 0;
}

sub diff_vcf {
    my ($self, $first_file, $second_file) = @_;
    my $first_md5  = qx(grep -vP '^##fileDate' $first_file | md5sum);
    my $second_md5 = qx(grep -vP '^##fileDate' $second_file | md5sum);
    return ($first_md5 eq $second_md5 ? 1 : 0);
}

sub diff_hq {
    my ($self, $first_file, $second_file) = @_;
    my $first_md5  = qx(grep -vP '^##fileDate' $first_file | grep -vP '^##startTime' | grep -vP '^##cmdline' | md5sum);
    my $second_md5 = qx(grep -vP '^##fileDate' $second_file | grep -vP '^##startTime' | grep -vP '^##cmdline' | md5sum);
    return ($first_md5 eq $second_md5 ? 1 : 0);
}

sub diff_vcf_gz {
    my ($self, $first_file, $second_file) = @_;
    my $first_md5  = qx(zcat $first_file | grep -vP '^##fileDate' | md5sum);
    my $second_md5 = qx(zcat $second_file | grep -vP '^##fileDate' | md5sum);
    return ($first_md5 eq $second_md5 ? 1 : 0);
}

# This method takes another build id and compares that build against this one. It gets
# a list of all the files in both builds and attempts to find pairs of corresponding
# files. The files/dirs listed in the files_ignored_by_diff and dirs_ignored_by_diff
# are ignored entirely, while files listed by regex_files_for_diff are retrieved
# using regex instead of a simple string eq comparison.
sub compare_output {
    my ($self, $other_build_id) = @_;
    my $build_id = $self->build_id;
    confess "Require build ID argument!" unless defined $other_build_id;
    my $other_build = Genome::Model::Build->get($other_build_id);
    confess "Could not get build $other_build_id!" unless $other_build;

    unless ($self->model_id eq $other_build->model_id) {
        confess "Builds $build_id and $other_build_id are not from the same model!";
    }
    unless ($self->class eq $other_build->class) {
        confess "Builds $build_id and $other_build_id are not the same type!";
    }

    # Create hashes for each build, keys are paths relative to build directory and
    # values are full file paths
    my (%file_paths, %other_file_paths);
    require Cwd;
    for my $file (@{$self->files_in_data_directory}) {
        my $abs_path = Cwd::abs_path($file);
        next unless $abs_path; # abs_path returns undef if a subdirectory of file does not exist
        $file_paths{$self->full_path_to_relative($file)} = $abs_path;
    }
    for my $other_file (@{$other_build->files_in_data_directory}) {
        $other_file_paths{$other_build->full_path_to_relative($other_file)} = Cwd::abs_path($other_file);
    }

    # Now cycle through files in this build's data directory and compare with
    # corresponding files in other build's dir
    my %diffs;
    FILE: for my $rel_path (sort keys %file_paths) {
        my $abs_path = delete $file_paths{$rel_path};
        warn "abs_path ($abs_path) does not exist\n" unless (-e $abs_path);
        my $dir = $self->full_path_to_relative(dirname($abs_path));

        next FILE if -d $abs_path;
        next FILE if $rel_path =~ /server_location.txt/;
        next FILE if grep { $dir =~ /$_/ } $self->dirs_ignored_by_diff;
        next FILE if grep { $rel_path =~ /$_/ } $self->files_ignored_by_diff;

        # Gotta check if this file matches any of the supplied regex patterns.
        # If so, find the one (and only one) file from the other build that
        # matches the same pattern
        my ($other_rel_path, $other_abs_path);
        REGEX: for my $regex ($self->regex_files_for_diff) {
            next REGEX unless $rel_path =~ /$regex/;

            #check for captures to narrow the search for the matching file
            if($regex =~ /\([^?].*?\)/) {
                my $modified_regex = $regex;
                for($rel_path =~ /$regex/) {
                    #replace captures with their found values from the original file
                    $modified_regex =~ s/\([^?].*?\)/$_/;
                }

                $regex = $modified_regex;
            }

            my @other_keys = grep { $_ =~ /$regex/ } sort keys %other_file_paths;
            if (@other_keys > 1) {
                $diffs{$rel_path} = "multiple files from $other_build_id matched file name pattern $regex\n" . join("\n", @other_keys);
                map { delete $other_file_paths{$_} } @other_keys;
                next FILE;
            }
            elsif (@other_keys < 1) {
                $diffs{$rel_path} = "no files from $other_build_id matched file name pattern $regex";
                next FILE;
            }
            else {
                $other_rel_path = shift @other_keys;
                $other_abs_path = delete $other_file_paths{$other_rel_path};
            }
        }

        # If file name doesn't match any regex, assume relative paths are the same
        unless (defined $other_rel_path and defined $other_abs_path) {
            $other_rel_path = $rel_path;
            $other_abs_path = delete $other_file_paths{$other_rel_path};
            unless (defined $other_abs_path) {
                $diffs{$rel_path} = "no file $rel_path from build $other_build_id";
                next FILE;
            }
        }

        # Check if the files end with a suffix that requires special handling. If not,
        # just do an md5sum on the files and compare
        my $diff_result = 0;
        my %matching_regex_for_custom_diff = $self->matching_regex_for_custom_diff($abs_path);
        if (keys %matching_regex_for_custom_diff > 1) {
            die "Path ($abs_path) matched multiple regex_for_custom_diff ('" . join("', '", keys %matching_regex_for_custom_diff) . "')!\n";
        }
        elsif (keys %matching_regex_for_custom_diff == 1) {
            my ($key) = keys %matching_regex_for_custom_diff;
            my $method = "diff_$key";
            # Check build class first
            if ($self->can($method)) {
                $diff_result = $self->$method($abs_path, $other_abs_path);
            }
            # then model class
            elsif ($self->model_class->can($method)) {
                $diff_result = $self->model_class->$method($abs_path, $other_abs_path);
            }
            # cannot find it..die
            else {
                die "Custom diff method ($method) not implemented in " . $self->class . " or ".$self->model_class."!\n";
            }
        }
        else {
            my $file_md5 = Genome::Sys->md5sum($abs_path);
            my $other_md5 = Genome::Sys->md5sum($other_abs_path);
            $diff_result = ($file_md5 eq $other_md5);
        }

        unless ($diff_result) {
            my $build_dir = $self->data_directory;
            my $other_build_dir = $other_build->data_directory;
            $diffs{$rel_path} = "files are not the same (diff -u {$build_dir,$other_build_dir}/$rel_path)";
        }
    }

    # Make sure the other build doesn't have any extra files
    for my $rel_path (sort keys %other_file_paths) {
        my $abs_path = delete $other_file_paths{$rel_path};
        warn "abs_path ($abs_path) does not exist\n" unless (-e $abs_path);
        my $dir = $self->full_path_to_relative(dirname($abs_path));
        next if -d $abs_path;
        next if grep { $dir =~ /$_/ } $self->dirs_ignored_by_diff;
        next if grep { $rel_path =~ /$_/ } $self->files_ignored_by_diff;
        $diffs{$rel_path} = "no file in build $build_id";
    }

    # Now compare metrics of both builds
    my %metric_diffs = $self->diff_metrics($other_build);
    @diffs{ keys %metric_diffs } = values %metric_diffs if %metric_diffs;

    return %diffs;
}

sub diff_metrics {
    my ($build1, $build2) = @_;

    my %diffs;
    my %metrics;
    map { $metrics{$_->name} = $_ } $build1->metrics;
    my %other_metrics;
    map { $other_metrics{$_->name} = $_ } $build2->metrics;

    METRIC: for my $metric_name (sort keys %metrics) {
        my $metric = $metrics{$metric_name};

        if ( grep { $metric_name =~ /$_/ } $build1->metrics_ignored_by_diff ) {
            delete $other_metrics{$metric_name} if exists $other_metrics{$metric_name};
            next METRIC;
        }

        my $other_metric = delete $other_metrics{$metric_name};
        unless ($other_metric) {
            $diffs{$metric_name} = "no build metric with name $metric_name found for build ".$build2->id;
            next METRIC;
        }

        my $metric_value = $metric->value;
        my $other_metric_value = $other_metric->value;
        unless ($metric_value eq $other_metric_value) {
            $diffs{$metric_name} = "metric $metric_name has value $metric_value for build ".$build1->id." and value " .
            "$other_metric_value for build ".$build2->id;
            next METRIC;
        }
    }

    # Catch any extra metrics that the other build has
    for my $other_metric_name (sort keys %other_metrics) {
        $diffs{$other_metric_name} = "no build metric with name $other_metric_name found for build ".$build1->id;
    }

    return %diffs;
}

sub input_differences_from_model {
    my $self = shift;

    my @build_inputs = $self->inputs;
    my @model_inputs = $self->model->inputs;

    #build a list of inputs to check against
    my %build_inputs;
    for my $build_input (@build_inputs) {
        $build_inputs{$build_input->name}{$build_input->value_class_name}{$build_input->value_id} = $build_input;
    }

    my @model_inputs_not_found;
    for my $model_input (@model_inputs) {
        my $build_input_found = delete($build_inputs{$model_input->name}{$model_input->value_class_name}{$model_input->value_id});

        unless ($build_input_found) {
            push @model_inputs_not_found, $model_input;
        }
    }

    my @build_inputs_not_found;
    for my $name (keys %build_inputs) {
        for my $value_class_name (keys %{ $build_inputs{$name} }) {
            for my $build_input_not_found (values %{ $build_inputs{$name}{$value_class_name} }) {
                my $value = $build_input_not_found->value;
                if($value and $value->isa('Genome::Model::Build') and $value->model and my ($model_input) = grep($_->value eq $value->model, @model_inputs_not_found) ) {
                    @model_inputs_not_found = grep($_ ne $model_input, @model_inputs_not_found);
                } else {
                    push @build_inputs_not_found, $build_input_not_found;
                }
            }
        }
    }

    return (\@model_inputs_not_found, \@build_inputs_not_found);
}

sub build_input_differences_from_model {
    return @{ ($_[0]->input_differences_from_model)[1] };
}

sub model_input_differences_from_model {
    return @{ ($_[0]->input_differences_from_model)[0] };
}

#a cheap convenience method for views
sub delta_model_input_differences_from_model {
    my $self = shift;

    my ($model_inputs, $build_inputs) = $self->input_differences_from_model;
    my @model_inputs_to_include;
    for my $model_input (@$model_inputs) {
        unless( grep{ $_->name eq $model_input->name } @$build_inputs ) {
            push @model_inputs_to_include, $model_input;
        }
    }
    return @model_inputs_to_include;
}

sub is_used_as_model_or_build_input {
    # Both models and builds have this method and as such it is currently duplicated.
    # We don't seem to have any place to put things that are common between Models and Builds.
    my $self = shift;

    my @model_inputs = Genome::Model::Input->get(
        value_id => $self->id,
        value_class_name => $self->class,
    );

    my @build_inputs = Genome::Model::Build::Input->get(
        value_id => $self->id,
        value_class_name => $self->class,
    );

    my @inputs = (@model_inputs, @build_inputs);

    return (scalar @inputs) ? 1 : 0;
}

sub child_workflow_instances {
    my $self = shift;
    return $self->_get_workflow_instance_children($self->newest_workflow_instance);
}

sub child_lsf_jobs {
    my $self = shift;
    my @workflow_instances = $self->child_workflow_instances;
    return unless @workflow_instances;
    my @dispatch_ids = grep {defined $_} map($_->current->dispatch_identifier, @workflow_instances);
    my @valid_ids = grep {$_ !~ /^P/} @dispatch_ids;
    return @valid_ids;
}

sub _get_workflow_instance_children {
    my $self = shift;
    my $parent = shift || return;
    return $parent, map($self->_get_workflow_instance_children($_), $parent->related_instances);
}

sub _preprocess_subclass_description {
    my ($class, $desc) = @_;
    my $build_subclass_name = $desc->{class_name};
    my $model_subclass_name = $build_subclass_name;
    $model_subclass_name =~ s/::Build//;
    my $model_subclass_meta = eval { $model_subclass_name->__meta__; };
    if ($model_subclass_meta) {
        my @model_properties = $model_subclass_meta->properties();
        my $has = $desc->{has};
        for my $p (@model_properties) {
            my $name = $p->property_name;
            if ($has->{$name}) {
                #warn "exists: $name on $build_subclass_name\n";
                next;
            }
            if (grep { $_->can($name) } (qw/Genome::Model Genome::Model::Build /, @{$desc->{is}}) ) {
                #warn "in parent: $name on $build_subclass_name\n";
                next;
            }
            
            my %data = %{$p};
            my $type = $data{data_type};
            if (!$type) {
                #warn "no type on $name for $model_subclass_name, cannot infer build properties\n";
                next;
            }
            for my $key (keys %data) {
                delete $data{$key} unless $key =~ /^is_/;
            }
            delete $data{is_specified_in_module_header};
            if ($type->isa("Genome::Model")) {
                $type =~ s/^Genome::Model/Genome::Model::Build/;
                $name =~ s/_model(?=($|s$))/_build/;
            }
            $data{data_type} = $type;
            $data{property_name} = $name;
            
            if (($p->can("is_input") and $p->is_input) or ($p->can("is_metric") and $p->is_metric)) {
                # code below will augment these with via/to/where
                #warn "build gets input/metric $name ($build_subclass_name)\n";
            }
            elsif ($data{via} and $data{via} eq 'last_complete_build') {
                # model properites which go through the last complete build exist directly on the build
                #%data = $data{meta_for_build_attribute}; 
                #warn "build gets direct $name ($build_subclass_name)\n";
            }
            #elsif (not $data{via} or ($data{via} ne 'inputs' and $data{via} ne 'metrics') ) {
            elsif ($data{via} and $data{via} =~ /^(subject)$/) {
                # all other model properties are accessible via the model
                #warn "build gets indirect $name ($build_subclass_name)\n";
                $data{via} = 'model';
                $data{to} = $name;
                $data{is_delegated} = 1;
                $data{is_calculated} = 1;
            }
            else {
                next;
            }

            $has->{$name} = \%data;
        }
    }
    
    my @names = keys %{ $desc->{has} };
    for my $prop_name (@names) {
        my $prop_desc = $desc->{has}{$prop_name};
       
        # skip old things for which the developer has explicitly set-up indirection
        next if $prop_desc->{id_by};
        next if $prop_desc->{via};
        next if $prop_desc->{reverse_as};
        next if $prop_desc->{implied_by};

        if ($prop_desc->{is_param} and $prop_desc->{is_input}) {
            die "class $class has is_param and is_input on the same property! $prop_name";
        }

        if (exists $prop_desc->{'is_param'} and $prop_desc->{'is_param'}) {
            $prop_desc->{'via'} = 'processing_profile',
            $prop_desc->{'to'} = $prop_name;
            $prop_desc->{'is_mutable'} = 0;
            $prop_desc->{'is_delegated'} = 1;
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
                    if ($model_subclass_name->can('_has_legacy_input_types') and $model_subclass_name->_has_legacy_input_types) {
                        # once old snapshots are gone we can backfill a second row per input with correct class names
                        # then once those are gone we can delete this logic and delete the old rows
                        push @where_class, value_class_name => 'UR::Value';
                    }
                    else {
                        push @where_class, value_class_name => $prop_class;
                    }
                }
            }

            $desc->{has}{$assoc} = {
                property_name => $assoc,
                implied_by => $prop_name,
                is => 'Genome::Model::Build::Input',
                reverse_as => 'build',
                where => [ name => $prop_name, @where_class ],
                is_mutable => $prop_desc->{is_mutable},
                is_optional => $prop_desc->{is_optional},
                is_many => 1, #$prop_desc->{is_many},
            };

            %$prop_desc = (%$prop_desc,
                via => $assoc,
                to => Genome::Model->_resolve_to_for_prop_desc($prop_desc),
            );
        }

        # Metrics
        if ( exists $prop_desc->{is_metric} and $prop_desc->{is_metric} ) {
            $prop_desc->{via} = 'metrics';
            $prop_desc->{where} = [ name => join(' ', split('_', $prop_name)) ];
            $prop_desc->{to} = 'value';
            $prop_desc->{is_delegated} = 1;
            $prop_desc->{is_mutable} = 1;
        }
    }

    my ($ext) = ($desc->{class_name} =~ /Genome::Model::Build::(.*)/);
    my $pp_subclass_name = 'Genome::ProcessingProfile::' . $ext;

    my $pp_data = $desc->{has}{processing_profile} = {};
    $pp_data->{data_type} = $pp_subclass_name;
    $pp_data->{via} = 'model';
    $pp_data->{to} = 'processing_profile';

    return $desc;
}

sub heartbeat {
    my $self = shift;
    my %options = @_;
    my %heartbeat = $self->_heartbeat;
    return ( $options{verbose} ? $heartbeat{message} : $heartbeat{is_ok} );
}

sub heartbeat_verbose {
    return $_[0]->heartbeat(verbose => 1);
}

sub _heartbeat {
    my $self = shift;

    my %heartbeat = (
        id => $self->id,
        status => $self->status,
        is_ok => 0,
    );
    if (grep { $heartbeat{status} eq $_ } ('Succeeded', 'Preserved')) {
        $heartbeat{is_ok} = 1;
        $heartbeat{message} = 'Build is succeeded. Stauts is '.$heartbeat{status}.'.';
        return %heartbeat;
    }

    unless (grep { $self->status eq $_ } ('Running', 'Scheduled')) {
        $heartbeat{message} = 'Build is not running/scheduled.';
        return %heartbeat;
    }

    my @wf_instances = ($self->newest_workflow_instance, $self->child_workflow_instances);
    my @wf_instance_execs = map { $_->current } @wf_instances;

    WF: for my $wf_instance_exec (@wf_instance_execs) {
        my $lsf_job_id = $wf_instance_exec->dispatch_identifier;
        my $wf_instance_exec_status = $wf_instance_exec->status;
        my $wf_instance_exec_id = $wf_instance_exec->execution_id;

        if (grep { $wf_instance_exec_status eq $_ } ('new', 'done')) {
            next WF;
        }

        # only certaion operation types would have LSF jobs and everything below is inspecting LSF status
        my $operation_type = $wf_instance_exec->operation_instance->operation->operation_type;
        unless ( grep { $operation_type->isa($_) } ('Workflow::OperationType::Command', 'Workflow::OperationType::Event') ) {
            next WF;
        }

        unless ($lsf_job_id) {
            $heartbeat{message} = "Workflow Instance Execution (ID: $wf_instance_exec_id) status ($wf_instance_exec_status) has no LSF job ID";
            last WF;
        }

        if ($lsf_job_id =~ /^P/) {
            next WF;
        }

        my $bjobs_output = qx(bjobs -l $lsf_job_id 2> /dev/null | tr '\\n' '\\0' | sed -r -e 's/\\x0\\s{21}//g' -e 's/\\x0/\\n\\n/g');
        chomp $bjobs_output;
        unless($bjobs_output) {
            $heartbeat{message} = "Expected bjobs (LSF ID: $lsf_job_id) output but received none.";
            last WF;
        }

        my $lsf_status = $self->status_from_bjobs_output($bjobs_output);
        if ($wf_instance_exec_status eq 'scheduled' && $lsf_status ne 'pend') {
            $heartbeat{message} = "Workflow Instance Execution (ID: $wf_instance_exec_id) status ($wf_instance_exec_status) does not match LSF status ($lsf_status)";
            last WF;
        }
        elsif ($wf_instance_exec_status eq 'scheduled' && $lsf_status eq 'pend') {
            next WF;
        }

        if ($wf_instance_exec_status eq 'running' && $lsf_status ne 'run') {
            $heartbeat{message} = "Workflow Instance Execution (ID: $wf_instance_exec_id) status ($wf_instance_exec_status) does not match LSF status ($lsf_status)";
            last WF;
        }

        if ($wf_instance_exec_status eq 'crashed' && ($lsf_status eq 'done' || $lsf_status eq 'exit')) {
            $heartbeat{message} = "Workflow Instance Execution (ID: $wf_instance_exec_id) crashed.";
            last WF;
        }

        if ($wf_instance_exec_status ne 'running' || $lsf_status ne 'run') {
            $heartbeat{message} = "Missing state ($wf_instance_exec_status/$lsf_status) condition, only running/run should reach this point";
            last WF;
        }

        my @pids = $self->pids_from_bjobs_output($bjobs_output);
        my $execution_host = $self->execution_host_from_bjobs_output($bjobs_output);
        unless ($execution_host) {
            $heartbeat{message} = 'Expected execution host.';
            last WF;
        }
        my $ps_cmd = "ssh $execution_host ps -o pid= -o stat= -p " . join(" -p ", @pids) . ' 2> /dev/null';
        my @ps_output = qx($ps_cmd);
        chomp(@ps_output);

        if (@ps_output != @pids) {
            $heartbeat{message} = 'Expected ps output for ' . @pids . ' PIDs (' . $execution_host . ': ' . join(', ', @pids) . ').';
            last WF;
        }

        for my $ps_output (@ps_output) {
            my ($stat) = $ps_output =~ /\d+\s+(.*)/;
            unless($stat =~ /^(R|S)/) {
                $heartbeat{message} = 'Expected PID to be in a R or S stat.';
                last WF;
            }
        }

        my $output_file = $wf_instance_exec->stdout;
        my $output_stat = stat($output_file);
        my $elapsed_mtime_output_file = time - $output_stat->mtime;
        my $error_file = $wf_instance_exec->stderr;
        my $error_stat = stat($error_file);
        my $elapsed_mtime_error_file = time - $error_stat->mtime;
        my $op_class = $wf_instance_exec->operation_instance->operation_type->command_class_name;
        my $hour = 3600;
        my $max_elapsed_time = ($op_class->can('max_elapsed_log_time')) ? $op_class->max_elapsed_log_time : 48 * $hour;
        if (($elapsed_mtime_output_file > $max_elapsed_time) && ($elapsed_mtime_error_file > $max_elapsed_time)) {
            my $elapsed_mtime_output_file_hours = int($elapsed_mtime_output_file/$hour);
            my $elapsed_mtime_error_file_hours = int($elapsed_mtime_error_file/$hour);
            my $max_elapsed_time_hours = int($max_elapsed_time/$hour);
            my $m = "Process is running BUT output and/or error file have not been modified in %d+ hours (%d hours, %d hours):\nOutput File: %s\nError File: %s";
            $heartbeat{message} = sprintf($m, $max_elapsed_time_hours, $elapsed_mtime_output_file_hours, $elapsed_mtime_error_file_hours, $output_file, $error_file);

            last WF;
        }
        $heartbeat{message} = 'OK. Seems to be running!';
        $heartbeat{is_ok} = 1;
    }

    return %heartbeat;
}

sub status_from_bjobs_output {
    my $self = shift;
    my $bjobs_output = shift;
    my ($status) = $bjobs_output =~ /Status <(.*?)>/;
    unless ($status) {
        die $self->error_message("Failed to parse status from bjobs output:\n$bjobs_output\n");
    }
    return lc($status);
}

sub pids_from_bjobs_output {
    my $self = shift;
    my $bjobs_output = shift;
    my ($pids) = $bjobs_output =~ /PIDs:([\d\s]+)/;
    unless ($pids) {
        die $self->error_message("Failed to parse PIDs from bjobs output:\n$bjobs_output\n");
    }
    my @pids = $pids =~ /(\d+)/;
    return @pids;
}

sub execution_host_from_bjobs_output {
    my $self = shift;
    my $bjobs_output = shift;
    my ($execution_host) = $bjobs_output =~ /Started on <(.*?)>/;
    unless ($execution_host) {
        if (my ($n_hosts, $hosts) = $bjobs_output =~ /Started on (\d+) Hosts\/Processors <(\S+)>/) {
            if ($hosts =~ /></) {
                $self->error_message("Not yet able to parse multiple execution hosts.");
            } else {
                ($execution_host) = $hosts =~ /$n_hosts\*(.*)/;
                $execution_host = $hosts unless $execution_host;
            }
        } else {
            die $self->error_message("Failed to parse execution host from bjobs output:\n$bjobs_output\n");
        }
    }
    my $ip_addr = qx(dig +short $execution_host); chomp $ip_addr;
    unless ($ip_addr) {
        $self->error_message("execution_host not recognized in DNS ($execution_host)");
        $execution_host = undef;
    }
    unless ($ip_addr =~ /^10\./) {
        $self->error_message("execution_host not recognized in DNS ($execution_host => $ip_addr)");
        $execution_host = undef;
    }
    return $execution_host;
}

sub is_current {
    my $self = shift;
    my $model = $self->model;

    my @build_inputs = $self->inputs;
    my @model_inputs = $model->inputs;
    unless ($model->_input_counts_are_ok(scalar(@model_inputs), scalar(@build_inputs))) {
        return;
    }

    my ($model_inputs_not_found, $build_inputs_not_found) = $self->input_differences_from_model;
    if (@$model_inputs_not_found || @$build_inputs_not_found) {
        unless ($model->_input_differences_are_ok($model_inputs_not_found, $build_inputs_not_found)) {
            return;
        }
    }

    my @from_builds = $self->from_builds;
    for my $from_build (@from_builds) {
        unless ($from_build->is_current) {
            return;
        }
    }

    return 1;
}

1;
