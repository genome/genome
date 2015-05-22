package Genome::InstrumentData::Composite::Workflow::Generator;

use strict;
use warnings;

use Genome;
use Sort::Naturally qw(nsort);
use Genome::Utility::Text qw(string_to_camel_case);

class Genome::InstrumentData::Composite::Workflow::Generator {
    is => 'Genome::InstrumentData::Composite::Workflow::Generator::Base',
};

my %VERSIONS = (
    'v1' => {
        picard_version => '1.29',
        samtools_version => 'r599',
    },
    'v2' => {
        picard_version => '1.46',
        samtools_version => 'r963',
    },
    'v3' => {
        picard_version => '1.42',
        samtools_version => 'r599',
    },
    'v4' => {
        picard_version => '1.85',
        samtools_version => 'r982',
    },
    'v5' => {
        picard_version => '1.113',
        samtools_version => 'r982',
    },
    'v6' => {
        picard_version => '1.113',
        samtools_version => '0.1.19',
        bedtools_version => '2.17.0',
    },
);

sub generate {
    my $class = shift;
    my $tree = shift;
    my $input_data = shift;
    my $merge_group = shift;

    my ($index_operations, $index_inputs) = Genome::InstrumentData::Composite::Workflow::Generator::AlignerIndex->generate($tree, $input_data);

    my @inputs;
    my ($object_workflows, $merge_operations, $refinement_operations, $refiners) = ({}, {}, {});

    my $objects_by_group = $class->_alignment_objects($tree, $input_data, $merge_group);
    for my $group (keys %$objects_by_group) {
        my $alignment_objects = $objects_by_group->{$group};

        my ($next_object_workflows, $next_object_inputs);
        my $alignment_generator = 'Genome::InstrumentData::Composite::Workflow::Generator::' . string_to_camel_case($tree->{action}->[0]->{type});
        ($next_object_workflows, $next_object_inputs) = $alignment_generator->generate( $tree, $input_data, $alignment_objects);
        @$object_workflows{keys %$next_object_workflows} = values %$next_object_workflows;
        push @inputs, @$next_object_inputs;

        my ($next_merge_operation, $next_merge_inputs) = Genome::InstrumentData::Composite::Workflow::Generator::Merge->generate($tree, $group, $alignment_objects);
        if($next_merge_operation) {
            $merge_operations->{$group} = $next_merge_operation;
            push @inputs, @$next_merge_inputs;
        }

        my ($next_refinement_operations, $next_refinement_inputs, $next_refiners) = Genome::InstrumentData::Composite::Workflow::Generator::Refine->generate($tree, $input_data, $group, $alignment_objects);
        if(@$next_refinement_operations) {
            $refinement_operations->{$group} = $next_refinement_operations;
            push @inputs, @$next_refinement_inputs;
            $refiners = $next_refiners;
        }
    }

    push @inputs, @$index_inputs;

    return $class->_generate_master_workflow($index_operations, $object_workflows, $merge_operations, $refinement_operations, \@inputs, $objects_by_group, $refiners, $tree->{api_version}, $input_data, $merge_group);
}

sub _alignment_objects {
    my $class = shift;
    my $tree = shift;
    my $input_data = shift;
    my $merge_group = shift;

    my @instrument_data = $class->_load_instrument_data($tree, $input_data);
    my @actions = @{$tree->{action}};
    my $read_aligner_name = $actions[0]->{name};

    my $output_by_group = {};

    my $should_segment = $read_aligner_name ne 'imported';
    for my $i (@instrument_data) {
        my $instrument_data_output = $output_by_group->{ $class->_merge_group_for_alignment_object($merge_group, $i) } ||= [];
        my @params = ($i, $class->_instrument_data_params($i, $input_data));
        my @segments;
        if ($should_segment && $i->isa('Genome::InstrumentData::Imported') && $i->can('get_segments')) {
            @segments = $i->get_segments();
        }

        if(@segments > 0) {
            for my $seg (@segments) {
                push @$instrument_data_output, [
                    @params,
                    instrument_data_segment_type => $seg->{segment_type},
                    instrument_data_segment_id => $seg->{segment_id},
                ];
            }
        } else {
            push @$instrument_data_output, \@params;
        }
    }

    return $output_by_group;
}

sub _load_instrument_data {
    my $class = shift;
    my $tree = shift;
    my $input_data = shift;

    my $instrument_data_param = $tree->{data};
    my $instrument_data = $input_data->{$instrument_data_param};
    unless($instrument_data and @$instrument_data) {
        die $class->error_message('No data found for ' . $instrument_data_param);
    }

    return @$instrument_data;
}

sub _merge_group_for_alignment_object {
    my $class = shift;
    my $merge_group = shift;
    my $alignment_object = shift;

    if($merge_group eq 'all') {
        return 'all';
    }

    if(ref $alignment_object eq 'ARRAY') {
        #accept either the output of _alignment_objects or a straight instrument_data
        $alignment_object = $alignment_object->[0];
    }

    my $group_obj = $alignment_object->$merge_group;
    unless($group_obj) {
        die $class->error_message('Could not determine ' . $merge_group . ' for data ' . $alignment_object->[0]->__display_name__);
    }

   return $group_obj->id;
}

sub _instrument_data_params {
    my $class = shift;
    my $instrument_data = shift;
    my $input_data = shift;

    return (instrument_data_id => $instrument_data->id, instrument_data_filter => $input_data->{'filter_' . $instrument_data->id});
}

sub _generate_master_workflow {
    my $class = shift;
    my $index_operations = shift;
    my $object_workflows = shift;
    my $merge_operations = shift;
    my $refinement_operations = shift;
    my $inputs = shift;
    my $objects_by_group = shift;
    my $refiners = shift;
    my $api_version = shift;
    my $input_data = shift;
    my $merge_group = shift;

    my ($input_properties_list, $output_properties_list) = $class->_inputs_and_outputs_for_master_workflow($index_operations, $object_workflows, $merge_operations, $refinement_operations, $refiners);

    my $master_workflow = Workflow::Model->create(
        name => 'Master Alignment Dispatcher',
        input_properties => $input_properties_list,
        optional_input_properties => $input_properties_list,
        output_properties => $output_properties_list,
    );

    my @block_operation_inputs = map { 'index_' . $_->{index} . '_result' } values %$index_operations;
    my $block_operation = Workflow::Operation->create(
        name => 'Wait for aligner indicies to be built',
        operation_type => Workflow::OperationType::Block->create(
            properties => [@block_operation_inputs,'force_fragment'],
        ),
    );
    $block_operation->workflow_model($master_workflow);
    $class->_add_link_to_workflow($master_workflow,
        left_workflow_operation_id => $master_workflow->get_input_connector->id,
        left_property => "m_force_fragment",
        right_workflow_operation_id => $block_operation->id,
        right_property => "force_fragment",
    );

    for my $index_operation (values %$index_operations){
        $class->_wire_index_operation_to_master_workflow($master_workflow, $block_operation, $index_operation);
    }

    for my $workflow (values %$object_workflows) {
        $class->_wire_object_workflow_to_master_workflow($master_workflow, $block_operation, $workflow);
    }

    if(%$merge_operations) {
        for my $merge_operation (values %$merge_operations) {
            $class->_wire_merge_operation_to_master_workflow($master_workflow, $block_operation, $merge_operation, $refinement_operations);
        }

        $class->_wire_object_workflows_to_merge_operations($master_workflow, $object_workflows, $merge_operations, $objects_by_group, $merge_group);
    }

    if(%$refinement_operations) {
        for my $group (keys %$refinement_operations) {
            my $refinement_operation = $refinement_operations->{$group};
            $class->_wire_refinement_operation_to_master_workflow($master_workflow, $block_operation, $refinement_operation, $refiners);
        }

        $class->_wire_merge_to_refinement_operations($master_workflow, $merge_operations, $refinement_operations);

        if (scalar(@$refiners) > 1) {
            $class->_wire_refinement_to_refinement_operations($master_workflow, $refinement_operations);
        }
    }

    #add the global inputs
    my $api_inputs = $class->inputs_for_api_version($api_version);
    my %inputs = (%$api_inputs, %$input_data);
    for my $input ($class->_general_workflow_input_properties) {
        push @$inputs,
            'm_' . $input => $inputs{$input};
    }

    return $master_workflow, $inputs;
}

sub inputs_for_api_version {
    my $class = shift;
    my $version = shift;

    unless(exists $VERSIONS{$version}) {
        die($class->error_message('Unknown API version: ' . $version));
    }

    return $VERSIONS{$version};
}

sub available_api_versions {
    return nsort keys %VERSIONS;
}

sub _inputs_and_outputs_for_master_workflow {
    my $class = shift;
    my $index_operations = shift;
    my $object_workflows = shift;
    my $merge_operations = shift;
    my $refinement_operations = shift;
    my $refiners = shift;

    my %input_properties;
    my %output_properties;

    for my $workflow (values %$object_workflows) {
        for my $prop (@{ $workflow->operation_type->input_properties }) {
            $input_properties{$prop}++;
        }
        for my $prop (@{ $workflow->operation_type->output_properties }) {
            $output_properties{$prop}++;
        }
    }

    if(%$merge_operations) {
        for my $prop ($class->_merge_workflow_input_properties) {
            $input_properties{$prop}++;
        }

        unless (%$refinement_operations) {
            for my $op (values %$merge_operations) {
                for my $prop (@{ $op->operation_type->output_properties }) {
                   $output_properties{join('_', $prop, $op->name)}++;
                }
            }
        }
    }

    if(%$refinement_operations) {
        for my $prop ($class->_refinement_workflow_input_properties($refiners)) {
            $input_properties{$prop}++;
        }

        my @last_refinement_operations = map { $_->[-1] } values %$refinement_operations;
        for my $op (@last_refinement_operations) {
            for my $prop (@{ $op->operation_type->output_properties }) {
               $output_properties{join('_', $prop, $op->name)}++;
            }
        }
    }

    #just need one copy if same parameter used in multiple subworkflows
    my @input_properties_list = map { 'm_' . $_ } keys %input_properties;
    my @output_properties_list = map { 'm_' . $_ } keys %output_properties;

    for my $index_operation (values %$index_operations) {
        push @input_properties_list, map { join('_', 'index', $index_operation->{index}, $_) } @{ $index_operation->{operation}->operation_type->input_properties };
    }

    return (\@input_properties_list, \@output_properties_list);
}

sub _wire_index_operation_to_master_workflow {
    my $class = shift;
    my $master_workflow = shift;
    my $block_operation = shift;
    my $operation = shift;

    $operation->{operation}->workflow_model($master_workflow);
        my $input_connector = $master_workflow->get_input_connector;
        my $output_connector = $master_workflow->get_output_connector;

        my $cname = $operation->{operation}->operation_type->command_class_name;
        my $cmeta = $cname->__meta__;
        for my $property (@{$operation->{operation}->operation_type->input_properties}){
            # This was written before the id_by properties were supported as inputs
            # so they were previously ignored.  We ignore them explicitly here
            # in this link autogen code because the associated ID property is
            # assigned directly.
            # In an updated version of the code we might skip things with implied_by
            # set to true, which would go the other way and capture the object,
            # while ignoring its supporting identity property.
            my $pmeta = $cmeta->property($property);
            if ($pmeta->id_by) {
                next;
            }

            if($property eq 'result_users') {
                #result users are the same for all steps in workflow
                $class->_add_link_to_workflow($master_workflow,
                    left_workflow_operation_id => $input_connector->id,
                    left_property => 'm_' . $property,
                    right_workflow_operation_id => $operation->{operation}->id,
                    right_property => $property,
                );
                next;
            }

            $class->_add_link_to_workflow($master_workflow,
                left_workflow_operation_id => $input_connector->id,
                left_property => join('_', "index", $operation->{index}, $property),
                right_workflow_operation_id => $operation->{operation}->id,
                right_property => $property,
            );
        }

        $class->_add_link_to_workflow($master_workflow,
            left_workflow_operation_id => $operation->{operation}->id,
            left_property => "result",
            right_workflow_operation_id => $block_operation->id,
            right_property => join("_", "index", $operation->{index}, "result"),
        );

    return 1;
}

sub _wire_object_workflow_to_master_workflow {
    my $class = shift;
    my $master_workflow = shift;
    my $block_operation = shift;
    my $workflow = shift;

    $workflow->workflow_model($master_workflow);

    #wire up the master to the inner workflows (just pass along the inputs and outputs)
    my $master_input_connector = $master_workflow->get_input_connector;
    my $workflow_input_connector = $workflow; #implicitly uses input connector
    for my $property (@{ $workflow->operation_type->input_properties }) {
        if($property eq 'force_fragment'){
            $class->_add_link_to_workflow($master_workflow,
                left_workflow_operation_id => $block_operation->id,
                left_property => $property,
                right_workflow_operation_id => $workflow_input_connector->id,
                right_property => $property,
            );
        }else {
            $class->_add_link_to_workflow($master_workflow,
                left_workflow_operation_id => $master_input_connector->id,
                left_property => 'm_' . $property,
                right_workflow_operation_id => $workflow_input_connector->id,
                right_property => $property,
            );
        }
    }

    my $master_output_connector = $master_workflow->get_output_connector;
    my $workflow_output_connector = $workflow; #implicitly uses output connector
    for my $property (@{ $workflow->operation_type->output_properties }) {
        $class->_add_link_to_workflow($master_workflow,
            left_workflow_operation_id => $workflow_output_connector->id,
            left_property => $property,
            right_workflow_operation_id => $master_output_connector->id,
            right_property => 'm_' . $property,
        );
    }

    return 1;
}

sub _wire_merge_operation_to_master_workflow {
    my $class = shift;
    my $master_workflow = shift;
    my $block_operation = shift; #unused by merge operations
    my $merge = shift;
    my $refinements = shift;

    my $master_input_connector = $master_workflow->get_input_connector;
    my $master_output_connector = $master_workflow->get_output_connector;

    $merge->workflow_model($master_workflow);

    for my $property ($class->_merge_workflow_input_properties) {
        $class->_add_link_to_workflow($master_workflow,
            left_workflow_operation_id => $master_input_connector->id,
            left_property => 'm_' . $property,
            right_workflow_operation_id => $merge->id,
            right_property => $property,
        );
    }

    unless (%$refinements) {
        for my $property (@{ $merge->operation_type->output_properties }) {
            $class->_add_link_to_workflow($master_workflow,
                left_workflow_operation_id => $merge->id,
                left_property => $property,
                right_workflow_operation_id => $master_output_connector->id,
                right_property => 'm_' . join('_', $property, $merge->name),
            );
        }
    }

    return 1;
}

sub _wire_refinement_operation_to_master_workflow {
    my $class = shift;
    my $master_workflow = shift;
    my $block_operation = shift; #unused by merge operations
    my $refinement_operation = shift;
    my $refiners = shift;

    my $master_input_connector = $master_workflow->get_input_connector;
    my $master_output_connector = $master_workflow->get_output_connector;

    for my $refinement (@$refinement_operation) {
        $refinement->workflow_model($master_workflow);
        my ($refiner) = grep { $refinement->name =~ /$_/ } @$refiners;

        for my $property ($class->_base_refinement_workflow_input_properties) {
            my $left_property = "m_" . $class->_construct_refiner_input_property($property, $refiner);
            my $right_property = $property;
            $class->_add_link_to_workflow($master_workflow,
                left_workflow_operation_id => $master_input_connector->id,
                left_property => $left_property,
                right_workflow_operation_id => $refinement->id,
                right_property => $right_property,
            );
        }

        $class->_add_link_to_workflow($master_workflow,
            left_workflow_operation_id => $master_input_connector->id,
            left_property => 'm_result_users',
            right_workflow_operation_id => $refinement->id,
            right_property => 'result_users',
        );
    }

    my $last_refinement = $refinement_operation->[-1];
    for my $property (@{ $last_refinement->operation_type->output_properties }) {
        $class->_add_link_to_workflow($master_workflow,
            left_workflow_operation_id => $last_refinement->id,
            left_property => $property,
            right_workflow_operation_id => $master_output_connector->id,
            right_property => 'm_' . join('_', $property, $last_refinement->name),
        );
    }

    return 1;
}

sub _wire_merge_to_refinement_operations {
    my $class = shift;
    my $master_workflow = shift;
    my $merge_operations = shift;
    my $refinement_operations = shift;

    for my $group_key (keys %$merge_operations ) {
        my $refinement_operation = $refinement_operations->{$group_key};
        $class->_add_link_to_workflow($master_workflow,
            left_workflow_operation_id => $merge_operations->{$group_key}->id,
            left_property => 'result_id',
            right_workflow_operation_id => $refinement_operation->[0]->id,
            right_property => 'input_result_id',
        );
    }

    return 1;
}

sub _wire_refinement_to_refinement_operations {
    my $class = shift;
    my $master_workflow = shift;
    my $refinement_operations = shift;

    for my $refinement_operation (values %$refinement_operations ) {
        next unless (scalar @$refinement_operation > 1);
        for my $i (1..$#$refinement_operation) {
            $class->_add_link_to_workflow($master_workflow,
                left_workflow_operation_id => $refinement_operation->[$i-1]->id,
                left_property => 'result_id',
                right_workflow_operation_id => $refinement_operation->[$i]->id,
                right_property => 'input_result_id',
            );
        }
    }

    return 1;
}


sub _wire_object_workflows_to_merge_operations {
    my $class = shift;
    my $master_workflow = shift;
    my $object_workflows = shift;
    my $merge_operations = shift;
    my $objects_by_group = shift;
    my $merge_group = shift;

    my %converge_inputs_for_group;
    for my $group (keys %$objects_by_group) {
        for my $o (@{ $objects_by_group->{$group} }) {
            my $align_wf = $object_workflows->{$o};

            $converge_inputs_for_group{$group} ||= [];
            push @{ $converge_inputs_for_group{$group} }, @{ $align_wf->operation_type->output_properties };
        }
    }

    my %converge_operation_for_group;
    for my $group (keys %$objects_by_group) {
        for my $o (@{ $objects_by_group->{$group} }) {
            my $align_wf = $object_workflows->{$o};
            my $merge_op = $merge_operations->{$group};

            unless(exists $converge_operation_for_group{$group}) {
                my $converge_operation = Workflow::Operation->create(
                    name => join('_', 'converge', $group),
                    operation_type => Workflow::OperationType::Converge->create(
                        input_properties => $converge_inputs_for_group{$group},
                        output_properties => ['alignment_result_ids'],
                    ),
                );

                $converge_operation->workflow_model($master_workflow);
                $class->_add_link_to_workflow($master_workflow,
                    left_workflow_operation_id => $converge_operation->id,
                    left_property => 'alignment_result_ids',
                    right_workflow_operation_id => $merge_op->id,
                    right_property => 'alignment_result_ids',
                );

                $converge_operation_for_group{$group} = $converge_operation;
            }

            for my $property (@{ $align_wf->operation_type->output_properties }) {
                $class->_add_link_to_workflow($master_workflow,
                    left_workflow_operation_id => $align_wf->id,
                    left_property => $property,
                    right_workflow_operation_id => $converge_operation_for_group{$group}->id,
                    right_property => $property,
                );
            }
        }
    }

    return 1;
}

sub _merge_workflow_input_properties {
    my $class = shift;

    return qw(
        merger_name
        merger_version
        merger_params
        duplication_handler_name
        duplication_handler_version
        duplication_handler_params
        bedtools_version
        samtools_version
        result_users
        );
}


sub _base_refinement_workflow_input_properties {
    my $class = shift;

    return qw(refiner_name refiner_version refiner_params refiner_known_sites_ids);
}

sub _refinement_workflow_input_properties {
    my $class = shift;
    my $refiners = shift;

    my @refiners = @$refiners;

    my @base_names = $class->_base_refinement_workflow_input_properties;

    my @input_properties;
    for my $refiner (@refiners) {
        for my $base_name (@base_names) {
            push @input_properties, $class->_construct_refiner_input_property($base_name, $refiner);
        }
    }

    return @input_properties;
}

1;
