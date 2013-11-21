
package Genome::ProcessingProfile::Staged;

use strict;
use warnings;
use Genome;
use Data::Dumper;

class Genome::ProcessingProfile::Staged {
    is => 'Genome::ProcessingProfile',
    is_abstract => 1,
    doc => 'processing profile subclass for workflows defined by events grouped into stages',
    has_param => [
        append_event_steps => {
            doc => 'Event classes to append to event_stage_job_classes, e.g. "alignment => Genome::Model::Event::Build::ReferenceAlignment::QC::CopyNumber".',
            is_optional => 1,
        },
    ],
};

sub stages {
    my $class = shift;
    $class = ref($class) if ref($class);
    die("Please implement stages in class '$class'");
}

sub classes_for_stage {
    my $self = shift;
    my $stage_name = shift;
    my $model = shift;
    my $classes_method_name = $stage_name .'_job_classes';
    my @classes = $self->$classes_method_name($model);
    push @classes, $self->_steps_to_append_from_processing_profile_for_stage($stage_name);
    return @classes;
}

sub objects_for_stage {
    my $self = shift;
    my $stage_name = shift;
    my $model = shift;
    my $objects_method_name = $stage_name .'_objects';
    return $self->$objects_method_name($model);
}

sub _steps_to_append_from_processing_profile_for_stage {
    my $self = shift;
    my $stage = shift;

    unless ($self->append_event_steps) {
        return;
    }

    my @stages = $self->stages;

    my $append_event_steps_value = $self->append_event_steps;
    my %append_event_steps;
    {
        no strict;
        %append_event_steps = eval("($append_event_steps_value)");
    };
    if (not %append_event_steps or $@) {
        my $err_msg = "Failed to interpret append_event_steps (" . $self->append_event_steps . ") as hash.";
        $err_msg .= "\nError: $@" if ($@);
        $err_msg .= "\n" . Dumper(\%append_event_steps);
        die $self->error_message($err_msg);
    }

    for my $key (keys %append_event_steps) {
        unless (grep { $_ =~ /$key/ } @stages) {
            die $self->error_message("Failed to find $key in stages (" . join(", ", @stages) . ").");
        }
    }

    if (exists $append_event_steps{$stage}) {
        my $steps = $append_event_steps{$stage};
        my @steps;
        if (ref $steps eq 'ARRAY') {
            @steps = @$steps;
        } elsif (my $ref = ref $steps) {
            die $self->error_message("Got REF ($ref), expects scalar string or ARRAY.");
        } else {
            @steps = ($steps);
        }

        my $pp_type = (split('::', $self->subclass_name))[2];
        for my $step (@steps) {
            unless ($step =~ /Genome\:\:Model\:\:Event\:\:Build\:\:$pp_type/) {
                die $self->error_message("Step ($step) does not appear to be appropriate for this processing profile type ($pp_type).");
            }

            my $class;
            {
                no strict;
                eval($class = $step->class);
            }
            unless($class) {
                die $self->error_message("Failed to load/find module '$step'.");
            }
            if($@) {
                die $self->error_message("Error ($@) while executing $step->class.");
            }
        }

        return @steps;
    } else {
        return;
    }
}

sub _resolve_workflow_for_build {
    # This is called by Genome::Model::Build::start()
    # Returns a Workflow::Operation
    # By default, builds this from stages(), but can be overridden for custom workflow.
    my $self = shift;
    my $build = shift;
    my $lsf_queue = shift; # TODO: the workflow shouldn't need this yet
    my $lsf_project = shift;

    if (!defined $lsf_queue || $lsf_queue eq '' || $lsf_queue eq 'inline') {
        $lsf_queue = $ENV{GENOME_LSF_QUEUE_BUILD_WORKER_ALT};
    }
    if (!defined $lsf_project || $lsf_project eq '') {
        $lsf_project = 'build' . $build->id;
    }

    my $events_by_stage = $self->_generate_events_for_build($build);

    my @workflow_stages;
    foreach my $stage_events ( @$events_by_stage ) {
        my $workflow_stage = $self->_workflow_for_stage( $build, $stage_events, $lsf_queue, $lsf_project )
            or next; # this is ok
        push @workflow_stages, $workflow_stage;
    }

    # ssmith ????  is this really ok?
    return unless @workflow_stages; # ok, may not have stages

    # FIXME check for errors here???
    my $workflow = $self->_merge_stage_workflows($build,@workflow_stages);

    return $workflow;
}

sub map_workflow_inputs {
    return (prior_result => 1);
}

sub _generate_events_for_build {
    my ($self, $build) = @_;

    my @stage_names = $self->stages($build);
    unless (@stage_names) {
        die ('No stages on processing profile in '  . (ref($self) || $self));
    }

    my @events_by_stage;
    for my $stage_name ( @stage_names ) {
        # FIXME why are we attempting to schedule stages that have no classes??
        my @events = $self->_generate_events_for_build_stage($build,$stage_name);
        unless ( @events ) {
            $self->warning_message(
                    'Stage '. $stage_name .' for build ('
                    . $build->build_id .") failed to schedule objects for classes:\n"
                    . join("\n",$self->classes_for_stage($stage_name, $build->model))
            );
            next;
        }
        push @events_by_stage, { name => $stage_name, events => \@events };
    }

    return $build->{events_by_stage} = \@events_by_stage;
}

sub _generate_events_for_build_stage {
    my ($self, $build, $stage_name) = @_;

    my @objects = $self->objects_for_stage($stage_name, $build->model);
    my @events;
    foreach my $object (@objects) {
        my $object_class;
        my $object_id;
        my $segment_identifier;
        if (ref($object)) {
            #if we get a segment passed in instead, extract the object & segment id
            if (ref($object) eq 'HASH' && exists $object->{segment}) {
                $segment_identifier = $object->{segment};
                $object = $object->{object};
            }
            $object_class = ref($object);
            $object_id = $object->id;
        } elsif ($object eq '1') {
            $object_class = 'single_instance';
        } else {
            $object_class = 'reference_sequence';
            $object_id = $object;
        }

        # Putting status message on build event because some tests expect it.
        #  Prolly can (re)move this to somewhere...
        if ($object_class->isa('Genome::InstrumentData')) {
            $build->debug_message('Scheduling jobs for '
                . $object_class . ' '
                . $object->full_name
                . ' (' . $object->id . ')'
            );
        } elsif ($object_class eq 'reference_sequence') {
            $build->debug_message('Scheduling jobs for reference sequence ' . $object_id);
        } elsif ($object_class eq 'single_instance') {
            $build->debug_message('Scheduling '. $object_class .' for stage '. $stage_name);
        } else {
            $build->debug_message('Scheduling for '. $object_class .' with id '. $object_id);
        }
        my @command_classes = $self->classes_for_stage($stage_name, $build->model);
        push @events, $self->_generate_events_for_object($build,$object,\@command_classes,$segment_identifier);
    }

    return @events;
}

sub _generate_events_for_object {
    my $self = shift;
    my $build = shift;
    my $object = shift;
    my $command_classes = shift;
    my $segment_identifier = shift;
    my $prior_event_id = shift;

    my @scheduled_commands;
    for my $command_class (@{$command_classes}) {
        if (ref($command_class) eq 'ARRAY') {
            push @scheduled_commands, $build->_schedule_command_classes_for_object($object,$command_class,$prior_event_id);
        } else {
            if ($command_class->can('command_subclassing_model_property')) {
                my $subclassing_model_property = $command_class->command_subclassing_model_property;
                unless ($build->model->$subclassing_model_property) {
                    # TODO: move into the creation of the processing profile
                    #$build->status_message("This processing profile doesNo value defined for $subclassing_model_property in the processing profile.  Skipping related processing...");
                    next;
                }
            }

            # Who did this????????????
            my $command;
            if ($command_class =~ /MergeAlignments|UpdateGenotype|FindVariations/) {
                if (ref($object)) {
                    unless ($object->isa('Genome::Model::RefSeq')) {
                        my $error_message = 'Expecting Genome::Model::RefSeq for EventWithRefSeq but got '. ref($object);
                        $build->error_message($error_message);
                        die;
                    }
                    $command = $command_class->create(
                        model_id => $build->model_id,
                        ref_seq_id => $object->ref_seq_id,
                    );
                } else {
                    $command = $command_class->create(
                        model_id => $build->model_id,
                        ref_seq_id => $object,
                    );
                }
            } elsif ($command_class =~ /ReferenceAlignment::AlignReads|ReferenceAlignment::BamQc|TrimReadSet|AssignReadSetToModel|AddReadSetToProject|FilterReadSet|RnaSeq::PrepareReads|DeNovoAssembly::ProcessInstrumentData/) {
                if ($object->isa('Genome::InstrumentData')) {
                    $command = $command_class->create(
                        instrument_data_id => $object->id,
                        model_id => $build->model_id,
                    );

                    if ($segment_identifier) {
                        $command->add_input(name=>'instrument_data_segment_type',
                                            value=>$segment_identifier->{segment_type});
                        $command->add_input(name=>'instrument_data_segment_id',
                                            value=>$segment_identifier->{segment_id})
                    }

                } else {
                    my $error_message = 'Expecting Genome::InstrumentData object but got '. ref($object);
                    $build->error_message($error_message);
                    die;
                }
            } elsif ($command_class->isa('Genome::Model::Event')) {
                $command = $command_class->create(
                    model_id => $build->model_id,
                    build_id => $build->id,
                );
            }
            unless ($command) {
                my $error_message = 'Problem creating subcommand for class '
                . ' for object class '. ref($object)
                . ' model id '. $build->model_id
                . ': '. $command_class->error_message();
                $build->error_message($error_message);
                die;
            }
            $command->build_id($build->build_id) unless defined $command->build_id;
            $command->prior_event_id($prior_event_id);
            $command->schedule;
            $prior_event_id = $command->id;
            push @scheduled_commands, $command;
            my $object_id;
            if (ref($object)) {
                $object_id = $object->id;
            } else {
                $object_id = $object;
            }
            $build->debug_message('Scheduled '. $command_class .' for '. $object_id
                .' event_id '. $command->genome_model_event_id ."\n");
        }
    }
    return @scheduled_commands;
}


sub _workflow_for_stage {
    my (
        $self,
        $build,
        $stage_from_build,
        $lsf_queue,          # TODO: this is passed from the build, but shouldn't be needed yet
        $lsf_project,
    ) = @_;

    my $stage_name = $stage_from_build->{name};
    my @events = @{$stage_from_build->{events}};
    unless (@events){
        $self->error_message('Failed to get events for stage '. $stage_name);
        return;
    }


    my $stage = Workflow::Model->create(
        name => $build->id . ' ' . $stage_name,
        input_properties => [
        'prior_result',
        ],
        output_properties => ['result']
    );
    my $input_connector = $stage->get_input_connector;
    my $output_connector = $stage->get_output_connector;

    my @ops_to_merge = ();
    my @first_events = grep { !defined($_->prior_event_id) } @events;
    for my $first_event ( @first_events ) {
        my $first_operation = $stage->add_operation(
            name => $first_event->command_name_brief .' '. $first_event->id,
            operation_type => Workflow::OperationType::Event->get(
                $first_event->id
            )
        );
        my $first_event_log_resource = $self->_resolve_log_resource($first_event);

        if (!$first_event_log_resource) {
            warn "no first_event_log_resource for first event " . $first_event->__display_name__;
            $first_event_log_resource = '';
        }

        if (!defined($first_event->bsub_rusage)) {
            print Data::Dumper::Dumper($first_event);
            Carp::confess();
        }

        $first_operation->operation_type->lsf_resource($first_event->bsub_rusage . $first_event_log_resource);
        $first_operation->operation_type->lsf_queue($lsf_queue);
        $first_operation->operation_type->lsf_project($lsf_project);

        $stage->add_link(
            left_operation => $input_connector,
            left_property => 'prior_result',
            right_operation => $first_operation,
            right_property => 'prior_result'
        );
        my $output_connector_linked = 0;
        my $sub;
        $sub = sub {
            my $prior_op = shift;
            my $prior_event = shift;
            my @events = $prior_event->next_events;
            if (@events) {
                foreach my $n_event (@events) {
                    my $n_operation = $stage->add_operation(
                        name => $n_event->command_name_brief .' '. $n_event->id,
                        operation_type => Workflow::OperationType::Event->get(
                            $n_event->id
                        )
                    );
                    my $n_event_log_resource = $self->_resolve_log_resource($n_event);
                    $n_operation->operation_type->lsf_resource($n_event->bsub_rusage . $n_event_log_resource);
                    $n_operation->operation_type->lsf_queue($lsf_queue);
                    $n_operation->operation_type->lsf_project($lsf_project);

                    $stage->add_link(
                        left_operation => $prior_op,
                        left_property => 'result',
                        right_operation => $n_operation,
                        right_property => 'prior_result'
                    );
                    $sub->($n_operation,$n_event);
                }
            } else {
                ## link the op's result to it.
                unless ($output_connector_linked) {
                    push @ops_to_merge, $prior_op;
                    $output_connector_linked = 1;
                }
            }
        };
        $sub->($first_operation,$first_event);
    }

    my $i = 1;
    my @input_names = map { 'result_' . $i++ } @ops_to_merge;
    my $converge = $stage->add_operation(
        name => 'merge results',
        operation_type => Workflow::OperationType::Converge->create(
            input_properties => \@input_names,
            output_properties => ['all_results','result']
        )
    );
    $i = 1;
    foreach my $op (@ops_to_merge) {
        $stage->add_link(
            left_operation => $op,
            left_property => 'result',
            right_operation => $converge,
            right_property => 'result_' . $i++
        );
    }
    $stage->add_link(
        left_operation => $converge,
        left_property => 'result',
        right_operation => $output_connector,
        right_property => 'result'
    );

    return $stage;
}

sub _merge_stage_workflows {
    my $self = shift;
    my $build = shift;
    my @workflows = @_;

    my %create_params = (
        name => $build->id . ' all stages',
        input_properties => [
            'prior_result'
        ],
        output_properties => [
            'result'
        ]
    );
    if ($ENV{WF_USE_FLOW}) {
        $create_params{'log_dir'} = sprintf("%s/logs/", $build->data_directory)
    }

    my $w = Workflow::Model->create(%create_params);

    my $last_op = $w->get_input_connector;
    my $last_op_prop = 'prior_result';
    foreach my $inner (@workflows) {
        $inner->workflow_model($w);

        $w->add_link(
            left_operation => $last_op,
            left_property => $last_op_prop,
            right_operation => $inner,
            right_property => 'prior_result'
        );

        $last_op = $inner;
        $last_op_prop = 'result';
    }

    $w->add_link(
        left_operation => $last_op,
        left_property => $last_op_prop,
        right_operation => $w->get_output_connector,
        right_property => 'result'
    );

    return $w;
}


1;

#$HeadURL$
#$Id$
