package Genome::WorkflowBuilder::DAG;

use strict;
use warnings;

use Genome;
use Params::Validate qw(:types);
use Set::Scalar qw();
use JSON;
use List::MoreUtils qw();
use Genome::Utility::Inputs qw(encode decode);
use Data::Dump qw(pp);
use File::Spec;
use Graph::Directed qw();
use Try::Tiny qw(try catch);


class Genome::WorkflowBuilder::DAG {
    is => 'Genome::WorkflowBuilder::Detail::Operation',

    has => [
        _optional_input_properties => {
            is => 'ARRAY',
            default => [],
        },
        operations => {
            is => 'ARRAY',
            default => [],
            doc => 'Genome::WorkflowBuilder::Detail::Operation objects',
        },

        links => {
            is => 'ARRAY',
            default => [],
            doc => 'Genome::WorkflowBuilder::Link objects',
        },

        _links => {
            is => 'HASH',
            default => {},
        },

        log_dir => {
            is => 'Text',
            is_optional => 1,
        },
    ],
    has_transient => {
        name_mapping_for_operation => {
            is => 'HASH',
            default => {},
        },
        outputs_for_operation => {
            is => 'HASH',
            default => {},
        }
    }
};

sub recursively_set_log_dir {
    my ($self, $log_dir) = Params::Validate::validate_pos(@_, 1,
        {type => SCALAR});

    $self->log_dir($log_dir);
    for my $op (@{$self->operations}) {
        if ($op->can('recursively_set_log_dir')) {
            $op->recursively_set_log_dir($log_dir);
        }
    }
    return;
}

sub parent_log_dir {
    my $class = shift;

    try {
        return Genome::Config::get('parent_workflow_log_directory');
    } catch {
        return;
    }
}

sub add_operation {
    my ($self, $op) = Params::Validate::validate_pos(@_, 1, {type => OBJECT});
    push @{$self->operations}, $op;

    my %constant_values = %{$op->constant_values};
    while (my ($constant_name, $value) = each %constant_values) {
        my $input_property = sprintf('%s.%s', $op->name, $constant_name);
        $self->connect_input(
            input_property => $input_property,
            destination => $op,
            destination_property => $constant_name,
        );
        $self->declare_constant($input_property => $value);
    }

    return $op;
}

sub add_link {
    my ($self, $link) = Params::Validate::validate_pos(@_, 1, 1);
    $self->_links->{$link->to_string} = $link;
    $self->links([values %{$self->_links}]);
    return $link;
}

# ------------------------------------------------------------------------------
# Public Methods
# ------------------------------------------------------------------------------

sub execute {
    my $self = shift;

    my %p = Params::Validate::validate(@_, {
        inputs => {type => HASHREF},
        polling_interval => {default => 120},
    });

    my $inputs = {%{$self->constant_values}, %{$p{inputs}}};

    my $backend = Genome::Config::get('workflow_builder_backend');
    if ($backend eq 'ptero') {
        return $self->_execute_with_ptero($inputs, $p{polling_interval});

    } elsif ($backend eq 'inline' or $backend eq 'simple') {
        return $self->execute_inline($inputs);

    } else {
        die sprintf("Unknown backend specified: %s", $backend);
    }
}

sub submit {
    my $self = shift;

    my %p = Params::Validate::validate(@_, {
        inputs => {type => HASHREF},
        process => {type => OBJECT},
    });

    my $inputs = {%{$self->constant_values}, %{$p{inputs}}};
    my $process = $p{process};

    my $backend = Genome::Config::get('workflow_builder_backend');
    if ($backend eq 'ptero') {
        my $wf_builder = $self->get_ptero_builder_for_process($process);

        my $wf_proxy = $wf_builder->submit(
            inputs => encode($inputs),
            submit_url => Genome::Config::get('ptero_workflow_submit_url'),
            name => $process->workflow_name
        );
        $self->status_message("Submitted workflow to petri service: %s",
            $wf_proxy->url);
        return $wf_proxy;
    } else {
        die sprintf("Only the ptero backend is supported " .
            "for 'submit' not: %s", $backend);
    }
}

sub _execute_with_ptero {
    my ($self, $inputs, $polling_interval) = @_;

    $self->remove_all_links_from_unspecified_inputs($inputs);
    my $used_inputs = $self->get_used_inputs($inputs);

    my $wf_builder = $self->get_ptero_builder($self->name);

    my $wf_proxy = $wf_builder->submit(
        inputs => encode($used_inputs),
        submit_url => Genome::Config::get('ptero_workflow_submit_url'),
    );
    $self->status_message("Waiting on PTero workflow (%s) to complete",
        $wf_proxy->url);
    Genome::Sys->disconnect_default_handles;
    $wf_proxy->wait(polling_interval => $polling_interval);

    if ($wf_proxy->has_succeeded) {
        if (!defined($wf_proxy->outputs)) {
            die sprintf('PTero workflow (%s) returned no results',
                $wf_proxy->url);
        } else {
            return decode($wf_proxy->outputs);
        }
    }
    else {
        die sprintf('PTero workflow (%s) did not succeed', $wf_proxy->url);
    }
}

sub get_ptero_builder {
    require Ptero::Builder::Workflow;

    my $self = shift;

    $self->validate;
    $self->remove_all_links_to_unused_input_properties();

    my $dag_method = Ptero::Builder::Workflow->new(name => $self->name);

    for my $operation (@{$self->operations}) {
        unless (defined $self->log_dir) {
            die sprintf("DAG (%s) has no log_dir set!", $self->name);
        }
        $dag_method->add_task($operation->get_ptero_builder_task($self->log_dir));
    }

    for my $link (@{$self->links}) {
        $link->validate;
        $dag_method->add_data_flow(
            source => $link->source_operation_name,
            source_property => $link->source_property,
            destination => $link->destination_operation_name,
            destination_property => $link->destination_property,
        );
    }

    return $dag_method;
}

sub get_ptero_builder_task {
    require Ptero::Builder::Detail::Workflow::Task;

    my $self = shift;

    $self->validate;

    my %params = (
        name => $self->name,
        methods => [
            $self->get_ptero_builder,
        ],
    );
    if (defined $self->parallel_by) {
        $params{parallel_by} = $self->parallel_by;
    }
    return Ptero::Builder::Detail::Workflow::Task->new(%params);
}


sub get_ptero_builder_for_process {
    require Ptero::Builder::Workflow;

    my ($self, $process) = Params::Validate::validate_pos(@_, 1, 1);
    my $process_id = $process->id;

    my $outer_dag = Ptero::Builder::Workflow->new(name => $process->workflow_name);

    my $set_to_running_task = $outer_dag->add_task(
        $self->_get_task_to_set_status($process_id, 'Running', 0)
    );
    $outer_dag->create_link(destination => $set_to_running_task);

    my $payload_task = $self->get_ptero_builder_task();
    $payload_task->methods([
        @{$payload_task->methods},
        $self->_get_method_to_set_status($process_id, 'Crashed', 1),
    ]);
    $outer_dag->add_task($payload_task);
    $outer_dag->create_link(
        source => $set_to_running_task,
        destination => $payload_task,
    );
    for my $input_property ($self->input_properties) {
        $outer_dag->add_data_flow(
            destination => $payload_task,
            source_property => $input_property,
            destination_property => $input_property,
        );
    }
    for my $output_property ($self->output_properties) {
        $outer_dag->add_data_flow(
            source => $payload_task,
            source_property => $output_property,
            destination_property => $output_property,
        );
    }

    my $set_to_succeeded_task = $outer_dag->add_task(
        $self->_get_task_to_set_status($process_id, 'Succeeded', 0)
    );
    $outer_dag->create_link(
        source => $payload_task,
        destination => $set_to_succeeded_task,
    );
    $outer_dag->create_link(
        source => $set_to_succeeded_task,
    );

    return $outer_dag;
}

sub _get_method_to_set_status {
    require Ptero::Builder::Job;

    my ($self, $process_id, $status, $exit_code) = Params::Validate::validate_pos(
        @_, 1, 1, 1, 1);

    return Ptero::Builder::Job->new(
        name => "set status $status",
        service_url => Genome::Config::get('ptero_shell_command_service_url'),
        parameters => {
            commandLine => [
                'genome', 'process', 'set-status',
                $process_id, $status,
                '--exit-code', $exit_code,
            ],
            environment => $self->_get_sanitized_env(),
            user => Genome::Sys->username,
            workingDirectory => Cwd::getcwd,
        },
    );
}

sub _get_task_to_set_status {
    require Ptero::Builder::Detail::Workflow::Task;

    my ($self, $process_id, $status, $exit_code) = Params::Validate::validate_pos(
        @_, 1, 1, 1, 1);

    my %params = (
        name => "set status $status",
        methods => [
            $self->_get_method_to_set_status($process_id, $status, $exit_code),
        ],
    );
    return Ptero::Builder::Detail::Workflow::Task->new(%params);
}

sub create_link {
    my $self = shift;
    $self->add_link(Genome::WorkflowBuilder::Link->create(@_));
    return;
}

sub connect_input {
    my $self = shift;
    my %args = Params::Validate::validate(@_, {
            input_property => { type => Params::Validate::SCALAR },
            destination => { type => Params::Validate::OBJECT },
            destination_property => { type => Params::Validate::SCALAR },
            is_optional => { type => Params::Validate::SCALAR, default => 0, },
    });

    $self->add_link(Genome::WorkflowBuilder::Link->create(
        source_property => $args{input_property},
        destination => $args{destination},
        destination_property => $args{destination_property},
    ));

    if ($args{is_optional}) {
        push @{$self->_optional_input_properties}, $args{input_property};
    }

    return;
}

sub connect_output {
    my $self = shift;
    my %args = Params::Validate::validate(@_, {
            source => { type => Params::Validate::OBJECT },
            source_property => { type => Params::Validate::SCALAR },
            output_property => { type => Params::Validate::SCALAR },
    });

    $self->add_link(Genome::WorkflowBuilder::Link->create(
        source => $args{source},
        source_property => $args{source_property},
        destination_property => $args{output_property},
    ));
    return;
}

sub operation_named {
    my ($self, $name) = @_;

    for my $op (@{$self->operations}) {
        if ($op->name eq $name) {
            return $op
        }
    }

    return;
}

sub is_input_property {
    my ($self, $property_name) = @_;

    return List::MoreUtils::any {$property_name eq $_} $self->input_properties;
}

sub is_optional_input_property {
    my ($self, $property_name) = @_;

    return List::MoreUtils::any {$property_name eq $_} $self->optional_input_properties;
}

sub is_output_property {
    my ($self, $property_name) = @_;

    return List::MoreUtils::any {$property_name eq $_} $self->output_properties;
}

sub is_many_property {
    my ($self, $property_name) = @_;
    # XXX There may not be an easy way to determine this.
    return;
}


# ------------------------------------------------------------------------------
# Inherited Methods
# ------------------------------------------------------------------------------

sub from_xml_element {
    my ($class, $element, $parent_log_dir) = @_;

    my $self = $class->create(
        name => $element->getAttribute('name'),
        log_dir => $element->getAttribute('logDir') || $parent_log_dir,
        parallel_by => $element->getAttribute('parallelBy'),
    );

    $self->_add_operations_from_xml_element($element);
    $self->_add_links_from_xml_element($element);
    $self->_add_optional_inputs_from_xml_element($element);

    return $self;
}

sub _add_optional_inputs_from_xml_element {
    my ($self, $element) = @_;

    my $ot_nodelist = $element->find('operationtype');
    for my $ot_node ($ot_nodelist->get_nodelist) {
        my $input_nodelist = $ot_node->find('inputproperty');
        for my $input_node ($input_nodelist->get_nodelist) {
            if ($input_node->getAttribute('isOptional')) {
                my $input_name = $input_node->textContent();
                push @{$self->_optional_input_properties}, $input_name;
            }
        }
    }
}

sub remove_all_links_from_unspecified_inputs {
    my ($self, $inputs) = @_;

    my @unspecified_inputs = $self->get_unspecified_optional_inputs($inputs);
    for my $input_property (@unspecified_inputs) {
        $self->remove_links_from_unspecified_input($input_property);
    }
    return;
}

sub remove_links_from_unspecified_input {
    my ($self, $input_property) = @_;

    if ($self->is_optional_input_property($input_property)) {
        my @kept_links;
        for my $link (@{$self->links}) {
            if ($link->external_input and
                $link->source_property eq $input_property) {
                my $destination = $link->destination;
                if ($destination) {
                    my $destination_property = $link->destination_property;
                    $destination->remove_links_from_unspecified_input(
                        $destination_property);
                } else {
                    push @kept_links, $link;
                }
            } else {
                push @kept_links, $link;
            }
        }
        $self->links([@kept_links]);
    } else {
        die sprintf("Cannot remove links from unspecified input (%s) " .
            "since it is not an optional input for DAG named (%s)",
            $input_property, $self->name);
    }
    return;
}

sub get_unspecified_optional_inputs {
    my ($self, $inputs) = @_;

    my @unspecified_inputs;
    for my $input_property ($self->optional_input_properties) {
        unless (exists $inputs->{$input_property}) {
            push @unspecified_inputs, $input_property;
        }
    }
    return @unspecified_inputs;
}

sub get_used_inputs {
    my ($self, $inputs) = @_;
    my %used_inputs;
    my @extra_inputs;
    for my $input_property (keys %$inputs) {
        if ($self->is_input_property($input_property)) {
            $used_inputs{$input_property} = $inputs->{$input_property};
        } else {
            if ($self->is_optional_input_property($input_property)) {
                $self->debug_message("Ignoring unused optional input named (%s)",
                    $input_property);
            } else {
                push @extra_inputs, $input_property;
            }
        }
    }

    if (@extra_inputs) {
        die sprintf("Extra inputs were specified to DAG->execute: %s",
            pp(\@extra_inputs));
    }

    return \%used_inputs;
}

sub remove_all_links_to_unused_input_properties {
    my $self = shift;
    $self->recurse_do('remove_links_to_unused_input_properties');
}

sub recurse_do {
    my ($self, $method_name) = @_;

    for my $operation (@{$self->operations}) {
        $operation->recurse_do($method_name);
    }
    $self->$method_name;
    return;
}

sub remove_links_to_unused_input_properties {
    my $self = shift;

    my $optionals = Set::Scalar->new($self->optional_input_properties);
    my @kept_links;
    for my $link (@{$self->links}) {
        if ($link->destination_is_unused_and_optional) {
            $self->debug_message("Removing link to DAG named (%s) that " .
                "targets unused optional input property named (%s)",
                $link->destination->name, $link->destination_property);
        } else {
            push @kept_links, $link;
        }
    }
    $self->links([@kept_links]);
    return;
}

sub get_xml_element {
    my $self = shift;

    my $element = $self->SUPER::get_xml_element($self);

    if (defined($self->log_dir)) {
        $element->setAttribute('logDir', $self->log_dir);
    }

    map {$element->addChild($_->get_xml_element)}
        sort {$a->name cmp $b->name} @{$self->operations};
    map {$element->addChild($_->get_xml_element)}
        sort {$a->sort_key cmp $b->sort_key} @{$self->links};

    return $element;
}

sub input_properties {
    my $self = shift;
    return sort $self->_property_names_from_links('external_input',
        'source_property');
}

sub optional_input_properties {
    my $self = shift;
    return sort @{$self->_optional_input_properties};
}

sub output_properties {
    my $self = shift;
    return sort $self->_property_names_from_links('external_output',
        'destination_property');
}

sub validate {
    my $self = shift;

    $self->SUPER::validate;

    $self->_validate_operation_names_are_unique;
    $self->_validate_linked_operation_ownership;
    $self->_validate_mandatory_inputs;
    $self->_validate_non_conflicting_inputs;

    for my $op (@{$self->operations}) {
        $op->validate;
    }

    for my $link (@{$self->links}) {
        $link->validate;
    }

    return;
}


# ------------------------------------------------------------------------------
# Private Methods
# ------------------------------------------------------------------------------

sub _add_operations_from_xml_element {
    my ($self, $element) = @_;

    my $nodelist = $element->find('operation');
    for my $node ($nodelist->get_nodelist) {
        my $op = Genome::WorkflowBuilder::Detail::Operation->from_xml_element(
            $node,
            $self->log_dir,
        );
        $self->add_operation($op);
    }
}

sub _add_links_from_xml_element {
    my ($self, $element) = @_;

    my $nodelist = $element->find('link');
    for my $node ($nodelist->get_nodelist) {
        my $source_op = $self->operation_named(
                $node->getAttribute('fromOperation'));
        my $destination_op = $self->operation_named(
                $node->getAttribute('toOperation'));

        my %link_params = (
            source_property => $node->getAttribute('fromProperty'),
            destination_property => $node->getAttribute('toProperty'),
        );
        if (defined($source_op)) {
            $link_params{source} = $source_op;
        }
        if (defined($destination_op)) {
            $link_params{destination} = $destination_op;
        }
        my $link = Genome::WorkflowBuilder::Link->create(%link_params);
        $self->add_link($link);
    }
}

sub _property_names_from_links {
    my ($self, $query_name, $property_holder) = @_;

    my $property_names = new Set::Scalar;

    for my $link (@{$self->links}) {
        if ($link->$query_name) {
            $property_names->insert($link->$property_holder);
        }
    }
    return @{$property_names};
}

sub _validate_operation_names_are_unique {
    my $self = shift;

    my $operation_names = new Set::Scalar;
    for my $op (@{$self->operations}) {
        if ($operation_names->contains($op->name)) {
            die sprintf("DAG '%s' contains multiple operations named '%s'",
                    $self->name, $op->name);
        }
        $operation_names->insert($op->name);
    }

    return;
}

sub _validate_linked_operation_ownership {
    my $self = shift;

    my %operations_hash;
    for my $op (@{$self->operations}) {$operations_hash{$op} = 1;}

    for my $link (@{$self->links}) {
        $self->_validate_operation_ownership($link->source, \%operations_hash);
        $self->_validate_operation_ownership($link->destination,
            \%operations_hash);
    }
    return;
}

sub _validate_operation_ownership {
    my ($self, $op, $operations_hash) = @_;

    if (defined($op)) {
        unless ($operations_hash->{$op}) {
            die sprintf("Unowned operation (%s) linked in DAG (%s)",
                    $op->name, $self->name,
            );
        }
    }
}

sub _validate_mandatory_inputs {
    my $self = shift;

    my $mandatory_inputs = $self->_get_mandatory_inputs;
    for my $link (@{$self->links}) {
        my $ei = $self->_encode_input($link->destination_operation_name,
            $link->destination_property);
        if ($mandatory_inputs->contains($ei)) {
            $mandatory_inputs->delete($ei);
        }
    }

    unless ($mandatory_inputs->is_empty) {
        die sprintf("%d mandatory input(s) missing in DAG: %s",
            $mandatory_inputs->size, $mandatory_inputs
        );
    }
}

sub _get_mandatory_inputs {
    my $self = shift;

    my $result = new Set::Scalar;

    for my $op (@{$self->operations}) {
        for my $prop_name ($op->input_properties) {
            $result->insert($self->_encode_input($op->name, $prop_name));
        }
    }

    return $result;
}

sub _encode_input {
    my ($self, $op_name, $property_name) = @_;
    my $js = JSON->new->allow_nonref;

    return $js->canonical->encode({
        operation_name => $op_name,
        property_name => $property_name,
    });
}

sub _validate_non_conflicting_inputs {
    my $self = shift;

    my $encoded_inputs = new Set::Scalar;
    for my $link (@{$self->links}) {
        my $ei = $self->_encode_input($link->destination_operation_name,
            $link->destination_property);
        if ($encoded_inputs->contains($ei)) {
            die sprintf("Conflicting input to '%s' on (%s) found.  " .
                "One link is from '%s' on (%s)",
                $link->destination_property, $link->destination_operation_name,
                $link->source_property, $link->source_operation_name
            );
        }
        $encoded_inputs->insert($ei);
    }
}

sub _execute_inline {
    my ($self, $inputs) = @_;

    $self->outputs_for_operation({});

    $self->_construct_name_mapping_for_operations();
    $self->_store_outputs_for_operation($inputs, "input connector");

    my @ordered_operations = $self->_get_ordered_operations();
    for my $operation (@ordered_operations) {
        my $inputs = $self->_get_inputs_for_operation($operation->name);

        my $outputs = $self->_execute_operation($operation, $inputs);

        $self->_store_outputs_for_operation($outputs, $operation->name);
    }

    return $self->_get_inputs_for_operation("output connector");
}

sub _execute_operation {
    my ($self, $operation, $inputs) = @_;

    if (defined($self->log_dir)) {
        Genome::Config::set_env('parent_workflow_log_directory', $self->log_dir);
    } else {
        $self->warning_message("No log directory set for DAG named (%s)", $self->name);
    }

    return $operation->execute_inline($inputs);
}

sub _construct_name_mapping_for_operations {
    my $self = shift;

    for my $link (@{$self->links}) {
        my $name_mapping = $self->name_mapping_for_operation;
        my $value = [$link->source_operation_name, $link->source_property];
        $name_mapping->{$link->destination_operation_name}->{$link->destination_property} = $value;
    }
    return;
}

sub _store_outputs_for_operation {
    my ($self, $outputs, $name) = @_;

    my $self_outputs = $self->outputs_for_operation;
    $self_outputs->{$name} = $outputs;
    return;
}

sub _get_ordered_operations {
    my $self = shift;

    my $g = Graph::Directed->new();
    for my $link (@{$self->links}) {
        $g->add_edge($link->source_operation_name, $link->destination_operation_name);
    }

    my @ordered_names;
    my @task_names = sort $g->successors('input connector');
    my $task_set = Set::Scalar->new(@task_names);
    $g->delete_vertex('input connector');
    while (scalar(@task_names) > 0) {
        my $count = 0;
        for my $name (@task_names) {
            if ($g->in_degree($name) == 0) {
                unless ($name eq 'output connector') {
                    push @ordered_names, $name;
                }
                splice(@task_names, $count, 1);

                my @new_successors = grep {!$task_set->contains($_)} $g->successors($name);
                push @task_names, sort @new_successors;
                $task_set->insert(@new_successors);

                $g->delete_vertex($name);
                last;
            }
            $count++;
        }
    }

    return map {$self->operation_named($_)} @ordered_names;
}

sub _get_inputs_for_operation {
    my ($self, $operation_name) = @_;

    my $inputs = {};
    my $name_mapping = $self->name_mapping_for_operation->{$operation_name};
    my $self_outputs = $self->outputs_for_operation;
    while (my ($input_name, $lookup_keys) = each(%{$name_mapping})) {
        $inputs->{$input_name} = $self_outputs->{$lookup_keys->[0]}->{$lookup_keys->[1]};
    }
    return $inputs;
}

1;
