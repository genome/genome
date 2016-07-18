package Genome::InstrumentData::Composite::Workflow::Generator::Align;

use strict;
use warnings;

use Genome;

class Genome::InstrumentData::Composite::Workflow::Generator::Align {
    is => 'Genome::InstrumentData::Composite::Workflow::Generator::Base',
};

sub generate {
    my $class = shift;
    my $master_workflow = shift;
    my $block_operation = shift;
    my $tree = shift;
    my $input_data = shift;
    my $alignment_objects = shift;

    my $workflows = {};
    my $inputs = [];

    for my $obj (@$alignment_objects) {
        my ($workflow, $input) = $class->_generate_workflow_for_instrument_data(
            $tree, $input_data, @$obj,
        );
        $workflows->{$obj} = $workflow;
        push @$inputs, @$input;

        $class->_wire_object_workflow_to_master_workflow($master_workflow, $block_operation, $workflow, $input);
    }


    return $workflows, $inputs;
}

sub _generate_workflow_for_instrument_data {
    my $class = shift;
    my $tree = shift;
    my $input_data = shift;
    my $instrument_data = shift;
    my %options = @_; #segment/read selection

    #First, get all the individual operations and figure out all the inputs/outputs we'll need
    my @operations = ();
    for my $subtree (@{ $tree->{action} } ) {
        push @operations, $class->_create_operations_for_alignment_tree($subtree, $instrument_data, %options);
    }

    #Next create the model, and add all the operations to it
    my $workflow_name = 'Alignment Dispatcher for ' . $instrument_data->id;
    if(exists $options{instrument_data_segment_id}) {
        $workflow_name .= ' (segment ' . $options{instrument_data_segment_id} . ')';
    }

    my $workflow = Genome::WorkflowBuilder::DAG->create(
        name => $workflow_name,
    );

    for my $op (@operations) {
        $workflow->add_operation($op);
    }

    #Last wire up the workflow
    my $inputs_for_links = $class->_generate_alignment_workflow_links($workflow, $tree, $input_data, $instrument_data, %options);

    push @$inputs_for_links, map { $class->_process_decorations($_, $workflow, $instrument_data, %options) } @{ $tree->{action} } ;

    return $workflow, $inputs_for_links;
}

sub _create_operations_for_alignment_tree {
    my $class = shift;
    my $tree = shift;
    my $instrument_data = shift;
    my %options = @_;

    my @operations = ();

    return @operations if exists $tree->{$class->_operation_key($instrument_data, %options)}; #generation already complete for this subtree

    push @operations, $class->_generate_operation($tree, $instrument_data, %options);
    if(exists $tree->{parent}) {
        push @operations, $class->_create_operations_for_alignment_tree($tree->{parent}, $instrument_data, %options);
    }

    return @operations;
}

sub _generate_alignment_workflow_links {
    my $class = shift;
    my $workflow = shift;
    my $tree = shift;
    my $input_data = shift;
    my $instrument_data = shift;
    my %options = @_;

    my $inputs = [];
    for my $subtree (@{ $tree->{action} }) {
        my $input = $class->_create_links_for_subtree($workflow, $subtree, $input_data, $instrument_data, %options);
        push @$inputs, @$input;

        #create link for final result of that subtree
        my $op = $subtree->{$class->_operation_key($instrument_data, %options)};
        $workflow->connect_output(
            source => $op,
            source_property => 'result_id',
            output_property => join('_', 'result_id', $op->name),
        );
    }

    return $inputs;
}

sub _create_links_for_subtree {
    my $class = shift;
    my $workflow = shift;
    my $tree = shift;
    my $input_data = shift;
    my $instrument_data = shift;
    my %options = @_;

    my $link_key = $class->_operation_key($instrument_data, %options) . '_links';
    return 1 if exists $tree->{$link_key};
    $tree->{$link_key} = 1;

    my $operation = $tree->{$class->_operation_key($instrument_data, %options)};
    die('No operation on tree!') unless $operation;

    my @properties = ('name', 'params', 'version');
    push @properties, 'reference_build_id'
        if ($tree->{type} eq 'align');
    push @properties, 'annotation_build_id'
        if (defined $tree->{annotation});

    my $inputs = [];
    for my $property (@properties) {
        my $source_property = join('_', $property, $operation->name);
        $workflow->connect_input(
            input_property => $source_property,
            destination => $operation,
            destination_property => $property,
        );

        my $value;

        if($property eq 'reference_build_id') {
            my $reference = $input_data->{$tree->{reference}};
            $value = $reference->id;
            unless($value) {
                $class->fatal_message('Found no reference for ' . $tree->{reference});
            }
        } elsif ($property eq 'annotation_build_id') {
            my $annotation_build = $input_data->{$tree->{annotation}};
            $value = $annotation_build->id;
            unless($value) {
                $class->fatal_message('Found no annotation build for ' . $tree->{annotation});
            }
        } else {
            $value = $tree->{$property};
        }

        push @$inputs, (
            'm_'.$source_property => $value,
        );
    }

    for my $property ($class->_general_workflow_input_properties) {
        $workflow->connect_input(
            input_property => $property,
            destination => $operation,
            destination_property => $property,
            is_optional => 1,
        );
    }

    if(exists $tree->{parent}) {
        my $parent_operation = $tree->{parent}{$class->_operation_key($instrument_data, %options)};
        $workflow->create_link(
            source => $parent_operation,
            source_property => 'output_result_id',
            destination => $operation,
            destination_property => 'alignment_result_input_id',
        );

        return $class->_create_links_for_subtree($workflow, $tree->{parent}, $instrument_data, %options);
    } else {
        #we're at the root--so create links for passing the segment info, etc.
        for my $property ($class->_instrument_data_workflow_input_properties($instrument_data, %options)) {
            my ($simple_property_name) = $property =~ m/^(.+?)__/;
            $workflow->connect_input(
                input_property => $property,
                destination => $operation,
                destination_property => $simple_property_name,
                is_optional => 1,
            );

            push @$inputs, (
                'm_'.$property => $options{$simple_property_name},
            );
        }
    }

    return $inputs;
}

sub _instrument_data_workflow_input_properties {
    my $class = shift;
    my $instrument_data = shift;
    my %options = @_;

    my @input_properties = map(
        join('__', $_, $instrument_data->id, (exists $options{instrument_data_segment_id} ? $options{instrument_data_segment_id} : ())),
        ('instrument_data_id', 'instrument_data_segment_id', 'instrument_data_segment_type', 'instrument_data_filter')
    );

    return @input_properties;
}

my $_generate_operation_tmpl;
sub _generate_operation {
    my $class = shift;
    my $action = shift;
    my $instrument_data = shift;
    my %options = @_;

    my $key = '_operation';
    if($instrument_data) {
        $key = $class->_operation_key($instrument_data, %options);
    }
    return $action->{$key} if exists $action->{$key}; #multiple leaves will point to the same parent, so stop if we've already worked on this subtree

    my $action_class_bases = {
        align => 'Genome::InstrumentData::Command::AlignReads',
        #filter => 'Genome::InstrumentData::Command::Filter', #TODO implement filter support
    };

    # TODO WF now has better support for customizing this, so this hack could be replaced
    # Use PerLaneTophat class to override default lsf_resource
    my $aligner_name = $action->{name};

    if ($aligner_name eq 'per-lane-tophat' or $aligner_name eq 'star') {
        my $subclass_name = 'Genome::InstrumentData::AlignmentResult'->_resolve_subclass_name_for_aligner_name($aligner_name);
        $action_class_bases->{align} = 'Genome::InstrumentData::Command::AlignReads::'.$subclass_name;
    }

    my $class_name = $action_class_bases->{$action->{type}};

    eval {
        $class_name->__meta__;
    };
    if(my $error = $@){
        $class->fatal_message("Could not create an instance of %s: %s", $class_name, $error);
    }

    unless ($_generate_operation_tmpl) {
        $_generate_operation_tmpl
            = UR::BoolExpr::Template->resolve('Genome::WorkflowBuilder::Command', 'id','name','command')->get_normalized_template_equivalent();
    }

    my $operation = Genome::WorkflowBuilder::Command->create(
        $_generate_operation_tmpl->get_rule_for_values(
            $class_name,
            UR::Object::Type->autogenerate_new_object_id_urinternal(),
            $class->get_unique_action_name($action),
        ),
    );

    my $aligner_class = join('::', 'Genome::InstrumentData::AlignmentResult', Genome::InstrumentData::AlignmentResult->_resolve_subclass_name_for_aligner_name($aligner_name));
    my $lsf_resource = $aligner_class->required_rusage(
        instrument_data => $instrument_data,
        aligner_params => $action->{params},
    );
    $operation->lsf_resource($lsf_resource);

    $action->{$key} = $operation;

    return $operation;
}

sub _operation_key {
    my $class = shift;
    my $instrument_data = shift;
    my %options = @_;

    my $key = '_operation_' . $instrument_data->id;
    if(exists $options{instrument_data_segment_id}) {
        $key .= '_' . $options{instrument_data_segment_id};
    }

    return $key;
}

sub get_unique_action_name {
    my $class = shift;
    my $action = shift;

    my $index = $class->next_counter_value;

    my $name = $action->{name};
    return join('_', $name, $index);
}

sub _process_decorations {
    my $class = shift;
    my $action = shift;
    my $workflow = shift;
    my $instrument_data = shift;
    my %options = @_;

    my @additional_inputs;

    if (exists $action->{decoration}) {
        my $operation_key = $class->_operation_key($instrument_data, %options);
        push @additional_inputs, Genome::InstrumentData::Composite::Decorator->decorate($action->{$operation_key}, $workflow, $action->{decoration});
    }

    if (exists $action->{parent}) {
        push @additional_inputs, $class->_process_decorations($action->{parent}, $instrument_data, %options);
    }

    return @additional_inputs;
}

sub _wire_object_workflow_to_master_workflow {
    my $class = shift;
    my $master_workflow = shift;
    my $block_operation = shift;
    my $workflow = shift;

    $master_workflow->add_operation($workflow);

    #wire up the master to the inner workflows (just pass along the inputs and outputs)
    for my $property ($workflow->input_properties) {
        if($property eq 'force_fragment'){
            $master_workflow->create_link(
                source => $block_operation,
                source_property => $property,
                destination => $workflow,
                destination_property => $property,
            );
        } else {
            $master_workflow->connect_input(
                input_property => 'm_' . $property,
                destination => $workflow,
                destination_property => $property,
                is_optional => 1,
            );
        }
    }

    for my $property ($workflow->output_properties) {
        $master_workflow->connect_output(
            source => $workflow,
            source_property => $property,
            output_property => 'm_' . $property,
        );
    }

    return 1;
}

1;
