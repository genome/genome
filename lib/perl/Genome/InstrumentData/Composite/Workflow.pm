package Genome::InstrumentData::Composite::Workflow;

use strict;
use warnings;

use Genome;

use Workflow;
use Workflow::Simple qw(run_workflow_lsf);

class Genome::InstrumentData::Composite::Workflow {
    is => 'Command',
    has => [
        inputs => {
            is => 'HASH',
            doc => 'Contains { name => [values] } for the terms in the alignment strategy',
        },
        strategy => {
            is => 'Text',
            doc => 'The alignment strategy to run',
        },
        merge_group => {
            is => 'Text',
            default_value => 'sample',
            valid_values => ['sample', 'all'],
            doc => 'How to group the instrument data when merging',
        },
    ],
    has_transient_optional => [
        _workflow => {
            is => 'Workflow::Model',
            doc => 'The underlying workflow to run the alignment/merge',
            is_output => 1,
        },
        _result_ids => {
            is_many => 1,
            doc => 'The alignments created/found as a result of running the workflow',
            is_output => 1,
        },
        log_directory => {
            is => 'Text',
            doc => 'Where to write the workflow logs',
        }
    ],
};

sub execute {
    my $self = shift;

    my $tree = $self->_process_strategy($self->strategy);
    die unless $tree;

    my ($workflow, $inputs) = $self->_generate_workflow($tree);
    $self->_workflow($workflow);

    if($self->log_directory) {
        $workflow->log_dir($self->log_directory);
    }

    my @errors = $workflow->validate;
    if (@errors) {
        $self->error_message(@errors);
        die 'Errors validating workflow';
    }

    my @results = $self->_run_workflow($workflow, $inputs);
    $self->_result_ids(\@results);

    return 1;
}

sub _process_strategy {
    my $self = shift;
    my $strategy_string = shift;

    $self->status_message('Analyzing strategy...');

    my $strategy = Genome::InstrumentData::Composite::Strategy->create(
        strategy => $strategy_string,
    );
    my $tree = $strategy->execute;

    unless( $tree ) {
        $self->error_message('Failed to analyze strategy: '.$self->strategy);
        return;
    }

    $self->status_message('Analyzing strategy...OK');
    return $tree;
}

sub _generate_workflow {
    my $self = shift;
    my $tree = shift;

    $self->status_message('Generating workflow...');

    unless(exists $tree->{data}) {
        die $self->error_message('No data specified in input');
    }

    my ($index_operations, $index_inputs) = $self->_generate_aligner_index_creation_operations($tree);
    my $instrument_data = $self->_load_instrument_data($tree);
    my @alignment_objects = $self->_alignment_objects($instrument_data, $tree);

    my ($object_workflows, $object_inputs) = $self->_generate_workflows_for_alignment_objects( $tree, \@alignment_objects);

    my ($merge_operations, $merge_inputs) = $self->_generate_merge_operations($tree, \@alignment_objects);

    my $inputs = [@$index_inputs, @$object_inputs, @$merge_inputs];

    return $self->_generate_master_workflow($index_operations, $object_workflows, $merge_operations, $inputs, \@alignment_objects, $tree->{api_version});
}

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
);

sub inputs_for_api_version {
    my $self = shift;
    my $version = shift;

    unless(exists $VERSIONS{$version}) {
        die($self->error_message('Unknown API version: ' . $version));
    }

    return $VERSIONS{$version};
}

sub _generate_aligner_index_creation_operations {
    my $self = shift;
    my $tree = shift;

    #find all of the aligners for which we need to check for indices
    my %aligners = ();

    my @nodes_to_visit = @{$tree->{action}};
    until(!@nodes_to_visit){
        my $current_node = shift @nodes_to_visit;
        if(exists $current_node->{parent}){
            push @nodes_to_visit, $current_node->{parent};
        }
        if($current_node->{type} eq 'align'){
            $aligners{$current_node->{name}}{$current_node->{version}}{$current_node->{reference}}{($current_node->{annotation} || '')}{$current_node->{params}}++;
        }
    }

    #build the workflow operations for each distinct index required
    my $workflow_operations = {};
    while(my ($aligner, $versions) = each %aligners) {
        my $alignment_result_class = 'Genome::InstrumentData::AlignmentResult::' . Genome::InstrumentData::AlignmentResult->_resolve_subclass_name_for_aligner_name($aligner);
        next unless $alignment_result_class && $alignment_result_class->can('prepare_reference_sequence_index');

        while(my ($version, $references) = each %$versions) {
            while(my ($reference, $annotations) = each %$references) {
                while(my ($annotation, $params_lists) = each %$annotations) {
                    if($alignment_result_class->aligner_params_required_for_index){
                        for my $params (keys %$params_lists) {
                            my $step_data = $self->_generate_aligner_index_creation_step($aligner, $version, $reference, $annotation, $params);
                            $workflow_operations->{$step_data->{operation}} = $step_data;
                        }
                    } else {
                        my $step_data = $self->_generate_aligner_index_creation_step($aligner, $version, $reference, $annotation);
                        $workflow_operations->{$step_data->{operation}} = $step_data;
                    }
                }
            }
        }
    }

    my @inputs;
    for my $property ($self->_index_workflow_input_properties) {
        for my $operation (values %$workflow_operations) {
            my $property_name = join('_', 'index', $operation->{index}, $property);
            push @inputs, $property_name => $operation->{$property};
        }
    }

    return ($workflow_operations, \@inputs);
}

sub _generate_aligner_index_creation_step {
    my $self = shift;
    my ($aligner, $version, $reference, $annotation, $params) = @_;

    my $index_num = $self->next_counter_value;
    my $operation = Workflow::Operation->create(
        name => "$aligner index #" . $index_num,
        operation_type => Workflow::OperationType::Command->create(
            command_class_name => 'Genome::Model::ReferenceSequence::Command::CreateAlignerIndex',
        ),
    );

    return {
        operation => $operation,
        index => $index_num,
        aligner_name => $aligner,
        aligner_version => $version,
        reference_sequence_build_id => $self->inputs->{$reference}->id,
        annotation_build_id => ($annotation? $self->inputs->{$annotation}->id : undef),
        aligner_params => (defined $params? $params : '')
 
    };
}

my $counter = 0;
sub next_counter_value {
    return ++$counter;
}

sub get_unique_action_name {
    my $self = shift;
    my $action = shift;

    my $index = $self->next_counter_value;

    my $name = $action->{name};
    return join('_', $name, $index);
}

sub _load_instrument_data {
    my $self = shift;
    my $tree = shift;

    my $instrument_data_param = $tree->{data};
    my $instrument_data = ($self->inputs)->{$instrument_data_param};
    unless($instrument_data and @$instrument_data) {
        die $self->error_message('No data found for ' . $instrument_data_param);
    }

    return $instrument_data;
}

sub _alignment_objects {
    my $self = shift;
    my $id = shift;
    my $tree = shift;
    my @instrument_data = @$id;
    my @actions = @{ $tree->{action}};
    my $read_aligner_name = $actions[0]->{name};

    my @instrument_data_output = map([$_, $self->_instrument_data_params($_)], grep {! $_->can('get_segments')} @instrument_data);
    my @segmentable_data = grep {$_->can('get_segments')} @instrument_data;

    for my $i (@segmentable_data) {
        my @segments = $i->get_segments();
        if (@segments > 1 && $read_aligner_name ne 'imported' && $i->isa('Genome::InstrumentData::Imported')) {
            for my $seg (@segments) {
                push @instrument_data_output, [
                    $i,
                    instrument_data_segment_type => $seg->{segment_type},
                    instrument_data_segment_id => $seg->{segment_id},
                    $self->_instrument_data_params($i),
                ];
            }
        } else {
            push @instrument_data_output, [$i, $self->_instrument_data_params($i)];
        }
    }

    return @instrument_data_output;
}

sub _instrument_data_params {
    my $self = shift;
    my $instrument_data = shift;

    return (instrument_data_id => $instrument_data->id, instrument_data_filter => ($self->inputs)->{'filter_' . $instrument_data->id});
}

sub _generate_workflows_for_alignment_objects {
    my $self = shift;
    my $tree = shift;
    my $alignment_objects = shift;

    my $workflows = {};
    my $inputs = [];

    for my $obj (@$alignment_objects) {
        my ($workflow, $input) = $self->_generate_workflow_for_instrument_data(
            $tree, @$obj,
        );
        $workflows->{$obj} = $workflow;
        push @$inputs, @$input;
    }

    return $workflows, $inputs;
}

sub _generate_workflow_for_instrument_data {
    my $self = shift;
    my $tree = shift;
    my $instrument_data = shift;
    my %options = @_; #segment/read selection

    #First, get all the individual operations and figure out all the inputs/outputs we'll need
    my @operations = ();
    for my $subtree (@{ $tree->{action} } ) {
        push @operations, $self->_create_operations_for_alignment_tree($subtree, $instrument_data, %options);
    }

    my @input_properties = (
        $self->_general_workflow_input_properties(),
        $self->_instrument_data_workflow_input_properties($instrument_data, %options),
        (map { $self->_input_properties_for_operation($_) } @operations)
    );

    my @output_properties;
    for my $leaf (@{ $tree->{action} }) {
        push @output_properties, join('_', 'result_id', $leaf->{$self->_operation_key($instrument_data, %options)}->name);
    }

    #Next create the model, and add all the operations to it
    my $workflow_name = 'Alignment Dispatcher for ' . $instrument_data->id;
    if(exists $options{instrument_data_segment_id}) {
        $workflow_name .= ' (segment ' . $options{instrument_data_segment_id} . ')';
    }

    my $workflow = Workflow::Model->create(
        name => $workflow_name,
        input_properties => \@input_properties,
        optional_input_properties => \@input_properties,
        output_properties => \@output_properties,
    );

    for my $op (@operations) {
        $op->workflow_model($workflow);
    }

    #Last wire up the workflow
    my $inputs_for_links = $self->_generate_alignment_workflow_links($workflow, $tree, $instrument_data, %options);

    return $workflow, $inputs_for_links;
}

sub _generate_alignment_workflow_links {
    my $self = shift;
    my $workflow = shift;
    my $tree = shift;
    my $instrument_data = shift;
    my %options = @_;

    my $inputs = [];
    for my $subtree (@{ $tree->{action} }) {
        my $input = $self->_create_links_for_subtree($workflow, $subtree, $instrument_data, %options);
        push @$inputs, @$input;

        #create link for final result of that subtree
        my $op = $subtree->{$self->_operation_key($instrument_data, %options)};
        $self->_add_link_to_workflow($workflow,
            left_workflow_operation_id => $op->id,
            left_property => 'result_id',
            right_workflow_operation_id => $workflow->get_output_connector->id,
            right_property => join('_', 'result_id', $op->name),
        );
    }

    return $inputs;
}

sub _create_links_for_subtree {
    my $self = shift;
    my $workflow = shift;
    my $tree = shift;
    my $instrument_data = shift;
    my %options = @_;

    my $link_key = $self->_operation_key($instrument_data, %options) . '_links';
    return 1 if exists $tree->{$link_key};
    $tree->{$link_key} = 1;

    my $operation = $tree->{$self->_operation_key($instrument_data, %options)};
    die('No operation on tree!') unless $operation;

    my @properties = ('name', 'params', 'version');
    push @properties, 'reference_build_id'
        if ($tree->{type} eq 'align');
    push @properties, 'annotation_build_id'
        if (defined $tree->{annotation});

    my $inputs = [];
    for my $property (@properties) {
        my $left_property = join('_', $property, $operation->name);
        $self->_add_link_to_workflow($workflow,
            left_workflow_operation_id => $workflow->get_input_connector->id,
            left_property => $left_property,
            right_workflow_operation_id => $operation->id,
            right_property => $property,
        );

        my $value;

        if($property eq 'reference_build_id') {
            my $reference = ($self->inputs)->{$tree->{reference}};
            $value = $reference->id;
            unless($value) {
                die $self->error_message('Found no reference for ' . $tree->{reference});
            }
        } elsif ($property eq 'annotation_build_id') {
            my $annotation_build = ($self->inputs)->{$tree->{annotation}};
            $value = $annotation_build->id;
            unless($value) {
                die $self->error_message('Found no annotation build for ' . $tree->{annotation});
            }
        } else {
            $value = $tree->{$property};
        }

        push @$inputs, (
            'm_'.$left_property => $value,
        );
    }

    for my $property ($self->_general_workflow_input_properties) {
        $self->_add_link_to_workflow($workflow,
            left_workflow_operation_id => $workflow->get_input_connector->id,
            left_property => $property,
            right_workflow_operation_id => $operation->id,
            right_property => $property,
        );
    }

    if(exists $tree->{parent}) {
        my $parent_operation = $tree->{parent}{$self->_operation_key($instrument_data, %options)};
        $self->_add_link_to_workflow($workflow,
            left_workflow_operation_id => $parent_operation->id,
            left_property => 'output_result_id',
            right_workflow_operation_id => $operation->id,
            right_property => 'alignment_result_input_id',
        );

        return $self->_create_links_for_subtree($workflow, $tree->{parent}, $instrument_data, %options);
    } else {
        #we're at the root--so create links for passing the segment info, etc.
        for my $property ($self->_instrument_data_workflow_input_properties($instrument_data, %options)) {
            my ($simple_property_name) = $property =~ m/^(.+?)__/;
            $self->_add_link_to_workflow($workflow,
                left_workflow_operation_id => $workflow->get_input_connector->id,
                left_property => $property,
                right_workflow_operation_id => $operation->id,
                right_property => $simple_property_name,
            );

            push @$inputs, (
                'm_'.$property => $options{$simple_property_name},
            );
        }
    }

    return $inputs;
}

my @_properties_for_add_link = ('id','workflow_model_id','left_workflow_operation_id','left_property',
                                'right_workflow_operation_id','right_property');
my @_properties_in_template_order;
sub _add_link_to_workflow {
    my($self,$workflow_obj, %params) = @_;

    unless ($self->{'_add_link_to_workflow_tmpl'}) {
        my $tmpl = UR::BoolExpr::Template->resolve('Workflow::Link', @_properties_for_add_link);
        unless ($tmpl) {
            Carp::croak('Cannot resolve BoolExpr template for adding Workflow::Link objects');
        }
        $tmpl = $tmpl->get_normalized_template_equivalent();
        for my $name ( @_properties_for_add_link ) {
            my $pos = $tmpl->value_position_for_property_name($name);
            if ($name eq 'workflow_model_id' ) {
                $self->{'_add_link_workflow_id_pos'} = $pos;
            } elsif ($name eq 'id') {
                $self->{'_add_link_id_pos'} = $pos;
            }
            $_properties_in_template_order[$pos] = $name;
        }
        $self->{'_add_link_to_workflow_tmpl'} = $tmpl;
    }
    my @values = @params{@_properties_in_template_order};
    $values[$self->{'_add_link_workflow_id_pos'}] = $workflow_obj->id;
    $values[$self->{'_add_link_id_pos'}] = UR::Object::Type->autogenerate_new_object_id_urinternal();
    my $rule = $self->{'_add_link_to_workflow_tmpl'}->get_rule_for_values(@values);
    return UR::Context->create_entity('Workflow::Link', $rule);
}

sub _instrument_data_workflow_input_properties {
    my $self = shift;
    my $instrument_data = shift;
    my %options = @_;

    my @input_properties = map(
        join('__', $_, $instrument_data->id, (exists $options{instrument_data_segment_id} ? $options{instrument_data_segment_id} : ())),
        ('instrument_data_id', 'instrument_data_segment_id', 'instrument_data_segment_type', 'instrument_data_filter')
    );

    return @input_properties;
}

sub _general_workflow_input_properties {
    my $self = shift;

    return qw(picard_version samtools_version force_fragment trimmer_name trimmer_version trimmer_params);
}

sub _merge_workflow_input_properties {
    my $self = shift;

    return qw(merger_name merger_version merger_params duplication_handler_name duplication_handler_version duplication_handler_params refiner_name refiner_version refiner_params refiner_known_sites_ids samtools_version);
}

sub _index_workflow_input_properties {
    my $self = shift;

    return qw(aligner_name aligner_params aligner_version reference_sequence_build_id annotation_build_id);
}

sub _create_operations_for_alignment_tree {
    my $self = shift;
    my $tree = shift;
    my $instrument_data = shift;
    my %options = @_;

    my @operations = ();

    return @operations if exists $tree->{$self->_operation_key($instrument_data, %options)}; #generation already complete for this subtree

    push @operations, $self->_generate_operation($tree, $instrument_data, %options);
    if(exists $tree->{parent}) {
        push @operations, $self->_create_operations_for_alignment_tree($tree->{parent}, $instrument_data);
    }

    return @operations;
}

my $_generate_operation_tmpl;
sub _generate_operation {
    my $self = shift;
    my $action = shift;
    my $instrument_data = shift;
    my %options = @_;

    my $key = '_operation';
    if($instrument_data) {
        $key = $self->_operation_key($instrument_data, %options);
    }
    return $action->{$key} if exists $action->{$key}; #multiple leaves will point to the same parent, so stop if we've already worked on this subtree

    my $action_class_bases = {
        align => 'Genome::InstrumentData::Command::AlignReads',
        #filter => 'Genome::InstrumentData::Command::Filter', #TODO implement filter support
    };

    # TODO WF now has better support for customizing this, so this hack could be replaced
    # Use PerLaneTophat class to override default lsf_resource
    if ($action->{name} eq 'per-lane-tophat') {
        $action_class_bases->{align} = 'Genome::InstrumentData::Command::AlignReads::PerLaneTophat';
    }

    my $class_name = $action_class_bases->{$action->{type}};

    eval {
        $class_name->__meta__;
    };
    if($@){
        $self->error_message("Could not create an instance of " . $class_name);
        die $self->error_message;
    }

    unless ($_generate_operation_tmpl) {
        # 'id' will still be first, since they're sorted alpha order
#ccc
        $_generate_operation_tmpl
            = UR::BoolExpr::Template->resolve('Workflow::Operation', 'id','name','workflow_operationtype_id')->get_normalized_template_equivalent();
    }

    my $operation = Workflow::Operation->create(
        $_generate_operation_tmpl->get_rule_for_values(UR::Object::Type->autogenerate_new_object_id_urinternal(),
                                                       $self->get_unique_action_name($action),
                                                       Workflow::OperationType::Command->get($class_name)->id),
    );

    $action->{$key} = $operation;

    return $operation;
}

sub _operation_key {
    my $self = shift;
    my $instrument_data = shift;
    my %options = @_;

    my $key = '_operation_' . $instrument_data->id;
    if(exists $options{instrument_data_segment_id}) {
        $key .= '_' . $options{instrument_data_segment_id};
    }

    return $key;
}

sub _input_properties_for_operation {
    my $self = shift;
    my $operation = shift;

    my @input_properties = ();
    my $op_name = $operation->name;
    my $op_type = $operation->operation_type;

    for my $prop ('reference_build_id', 'annotation_build_id') {
        if(grep $_ eq $prop, @{ $op_type->input_properties }) {
            push @input_properties, join('_', $prop, $op_name);
        }
    }

    push @input_properties, join('_', 'name', $op_name), join('_', 'params', $op_name), join('_', 'version', $op_name);

    return @input_properties;
}

sub _generate_merge_operations {
    my $self = shift;
    my $tree = shift;
    my $alignment_objects = shift;

    my %merge_operations;
    my @inputs;

    if(exists $tree->{then}) {
        my $merge_tree = $tree->{then};
        push @inputs, (
            m_merger_name => $merge_tree->{name},
            m_merger_version => $merge_tree->{version},
            m_merger_params => $merge_tree->{params},
        );

        my $next_op = $merge_tree;
        while(exists $next_op->{then}) {
            $next_op = $next_op->{then};

            if($next_op->{type} eq 'deduplicate') {
                push @inputs, (
                    m_duplication_handler_name => $next_op->{name},
                    m_duplication_handler_params => $next_op->{params},
                    m_duplication_handler_version => $next_op->{version}
                );
            } elsif($next_op->{type} eq 'refine') {
                push @inputs, (
                    m_refiner_name => $next_op->{name},
                    m_refiner_params => $next_op->{params},
                    m_refiner_version => $next_op->{version},
                    m_refiner_known_sites_ids => [ map { $_->id } @{($self->inputs)->{$next_op->{known_sites}}} ],
                );
            }
        }

        if($self->merge_group eq 'all') {
            #this case is a simplification for efficiency
            $merge_operations{all} = $self->_generate_merge_operation($merge_tree, 'all');
        } else {
            my %groups;
            for my $o (@$alignment_objects) {
                $groups{ $self->_merge_group_for_alignment_object($o) }++;
            }

            %merge_operations = map (($_ => $self->_generate_merge_operation($merge_tree, $_)), sort keys %groups);
        }
    }

    return (\%merge_operations, \@inputs);
}

sub _merge_group_for_alignment_object {
    my $self = shift;
    my $alignment_object = shift;

    if($self->merge_group eq 'all') {
        return 'all';
    }

    if(ref $alignment_object eq 'ARRAY') {
        #accept either the output of _alignment_objects or a straight instrument_data
        $alignment_object = $alignment_object->[0];
    }

    my $grouping = $self->merge_group;
    my $group_obj = $alignment_object->$grouping;
    unless($group_obj) {
        die $self->error_message('Could not determine ' . $grouping . ' for data ' . $alignment_object->[0]->__display_name__);
    }

   return $group_obj->id;
}

my $_mergealignments_command_id;
my $_generate_merge_operation_tmpl;
sub _generate_merge_operation {
    my $self = shift;
    my $merge_tree = shift;
    my $grouping = shift;

    unless ($_generate_merge_operation_tmpl) {
        $_generate_merge_operation_tmpl = UR::BoolExpr::Template->resolve('Workflow::Operation', 'id','name','workflow_operationtype_id')->get_normalized_template_equivalent();
        $_mergealignments_command_id = Workflow::OperationType::Command->get('Genome::InstrumentData::Command::MergeAlignments')->id;
    }

    my $operation = Workflow::Operation->create(
        $_generate_merge_operation_tmpl->get_rule_for_values(
                UR::Object::Type->autogenerate_new_object_id_urinternal(),
                'merge_' . $grouping,
                $_mergealignments_command_id
        )

        #name => 'merge_' . $grouping,
        #operation_type => $_mergealignments_command_id
    );

    return $operation;
}

sub _generate_master_workflow {
    my $self = shift;
    my $index_operations = shift;
    my $object_workflows = shift;
    my $merge_operations = shift;
    my $inputs = shift;
    my $alignment_objects = shift;
    my $api_version = shift;

    my ($input_properties_list, $output_properties_list) = $self->_inputs_and_outputs_for_master_workflow($index_operations, $object_workflows, $merge_operations);

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
    $self->_add_link_to_workflow($master_workflow,
        left_workflow_operation_id => $master_workflow->get_input_connector->id,
        left_property => "m_force_fragment",
        right_workflow_operation_id => $block_operation->id,
        right_property => "force_fragment",
    );

    for my $index_operation (values %$index_operations){
        $self->_wire_index_operation_to_master_workflow($master_workflow, $block_operation, $index_operation);
    }

    for my $workflow (values %$object_workflows) {
        $self->_wire_object_workflow_to_master_workflow($master_workflow, $block_operation, $workflow);
    }

    if(%$merge_operations) {
        for my $merge_operation (values %$merge_operations) {
            $self->_wire_merge_operation_to_master_workflow($master_workflow, $block_operation, $merge_operation);
        }

        $self->_wire_object_workflows_to_merge_operations($master_workflow, $object_workflows, $merge_operations, $alignment_objects);
    }

    #add the global inputs
    my $api_inputs = $self->inputs_for_api_version($api_version);
    my %inputs = (%$api_inputs, %{ $self->inputs });
    for my $input ($self->_general_workflow_input_properties) {
        push @$inputs,
            'm_' . $input => $inputs{$input};
    }

    return $master_workflow, $inputs;
}

sub _inputs_and_outputs_for_master_workflow {
    my $self = shift;
    my $index_operations = shift;
    my $object_workflows = shift;
    my $merge_operations = shift;

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
        for my $prop ($self->_merge_workflow_input_properties) {
            $input_properties{$prop}++;
        }

        for my $op (values %$merge_operations) {
            for my $prop (@{ $op->operation_type->output_properties }) {
               $output_properties{join('_', $prop, $op->name)}++;
            }
        }
    }

    #just need one copy if same parameter used in multiple subworkflows
    my @input_properties_list = map { 'm_' . $_ } keys %input_properties;
    my @output_properties_list = map { 'm_' . $_ } keys %output_properties;

    for my $index_operation (values %$index_operations) {
        for my $property ($self->_index_workflow_input_properties) {
            push @input_properties_list, join('_', 'index', $index_operation->{index}, $property);
        }
    }

    return (\@input_properties_list, \@output_properties_list);
}

sub _wire_index_operation_to_master_workflow {
    my $self = shift;
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
            $self->_add_link_to_workflow($master_workflow,
                left_workflow_operation_id => $input_connector->id,
                left_property => join('_', "index", $operation->{index}, $property),
                right_workflow_operation_id => $operation->{operation}->id,
                right_property => $property,
            );
        }

        $self->_add_link_to_workflow($master_workflow,
            left_workflow_operation_id => $operation->{operation}->id,
            left_property => "result",
            right_workflow_operation_id => $block_operation->id,
            right_property => join("_", "index", $operation->{index}, "result"),
        );

    return 1;
}

sub _wire_object_workflow_to_master_workflow {
    my $self = shift;
    my $master_workflow = shift;
    my $block_operation = shift;
    my $workflow = shift;

    $workflow->workflow_model($master_workflow);

    #wire up the master to the inner workflows (just pass along the inputs and outputs)
    my $master_input_connector = $master_workflow->get_input_connector;
    my $workflow_input_connector = $workflow; #implicitly uses input connector
    for my $property (@{ $workflow->operation_type->input_properties }) {
        if($property eq 'force_fragment'){
            $self->_add_link_to_workflow($master_workflow,
                left_workflow_operation_id => $block_operation->id,
                left_property => $property,
                right_workflow_operation_id => $workflow_input_connector->id,
                right_property => $property,
            );
        }else {
            $self->_add_link_to_workflow($master_workflow,
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
        $self->_add_link_to_workflow($master_workflow,
            left_workflow_operation_id => $workflow_output_connector->id,
            left_property => $property,
            right_workflow_operation_id => $master_output_connector->id,
            right_property => 'm_' . $property,
        );
    }

    return 1;
}

sub _wire_merge_operation_to_master_workflow {
    my $self = shift;
    my $master_workflow = shift;
    my $block_operation = shift; #unused by merge operations
    my $merge = shift;

    my $master_input_connector = $master_workflow->get_input_connector;
    my $master_output_connector = $master_workflow->get_output_connector;

    $merge->workflow_model($master_workflow);

    for my $property ($self->_merge_workflow_input_properties) {
        $self->_add_link_to_workflow($master_workflow,
            left_workflow_operation_id => $master_input_connector->id,
            left_property => 'm_' . $property,
            right_workflow_operation_id => $merge->id,
            right_property => $property,
        );
    }

    for my $property (@{ $merge->operation_type->output_properties }) {
        $self->_add_link_to_workflow($master_workflow,
            left_workflow_operation_id => $merge->id,
            left_property => $property,
            right_workflow_operation_id => $master_output_connector->id,
            right_property => 'm_' . join('_', $property, $merge->name),
        );
    }

    return 1;
}

sub _wire_object_workflows_to_merge_operations {
    my $self = shift;
    my $master_workflow = shift;
    my $object_workflows = shift;
    my $merge_operations = shift;
    my $alignment_objects = shift;

    my %converge_inputs_for_group;
    for my $o (@$alignment_objects) {
        my $align_wf = $object_workflows->{$o};
        my $group = $self->_merge_group_for_alignment_object($o);

        $converge_inputs_for_group{$group} ||= [];

        push @{ $converge_inputs_for_group{$group} }, @{ $align_wf->operation_type->output_properties };
    }

    my %converge_operation_for_group;
    for my $o (@$alignment_objects) {
        my $align_wf = $object_workflows->{$o};
        my $group = $self->_merge_group_for_alignment_object($o);
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
            $self->_add_link_to_workflow($master_workflow,
                left_workflow_operation_id => $converge_operation->id,
                left_property => 'alignment_result_ids',
                right_workflow_operation_id => $merge_op->id,
                right_property => 'alignment_result_ids',
            );

            $converge_operation_for_group{$group} = $converge_operation;
        }

        for my $property (@{ $align_wf->operation_type->output_properties }) {
            $self->_add_link_to_workflow($master_workflow,
                left_workflow_operation_id => $align_wf->id,
                left_property => $property,
                right_workflow_operation_id => $converge_operation_for_group{$group}->id,
                right_property => $property,
            );
        }
    }


    return 1;
}

sub _run_workflow {
    my $self = shift;
    my $workflow = shift;
    my $inputs = shift;

    $self->status_message('Running workflow...');
    my $result = Workflow::Simple::run_workflow_lsf( $workflow, @$inputs);

    unless($result){
        $self->error_message( join("\n", map($_->name . ': ' . $_->error, @Workflow::Simple::ERROR)) );
        die $self->error_message("Workflow did not return correctly.");
    }

    my @result_ids;
    for my $key (keys %$result) {
        if($key =~ /result_id/) {
            push @result_ids, $result->{$key};
        }
    }

    $self->status_message('Produced results: ' . join(', ', @result_ids));
    $self->status_message('Workflow complete.');

    return @result_ids;
}
1;
