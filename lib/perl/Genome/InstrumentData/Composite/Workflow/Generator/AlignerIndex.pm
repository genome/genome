package Genome::InstrumentData::Composite::Workflow::Generator::AlignerIndex;

use strict;
use warnings;

use Genome;

class Genome::InstrumentData::Composite::Workflow::Generator::AlignerIndex {
    is => 'Genome::InstrumentData::Composite::Workflow::Generator::Base',
};

sub generate {
    my $class = shift;
    my $master_workflow = shift;
    my $tree = shift;
    my $inputs = shift;

    #find all of the aligners for which we need to check for indices
    my %aligners = ();

    my @nodes_to_visit = @{$tree->{action}};
    until(!@nodes_to_visit){
        my $current_node = shift @nodes_to_visit;
        if(exists $current_node->{parent}){
            push @nodes_to_visit, $current_node->{parent};
        }
        if($current_node->{type} eq 'align' || $current_node->{type} eq 'align_and_merge'){
            $aligners{$current_node->{name}}{$current_node->{version}}{$current_node->{reference}}{($current_node->{annotation} || '')}{$current_node->{params}}++;
        }
    }

    #build the workflow operations for each distinct index required
    my $workflow_operations = {};
    while(my ($aligner, $versions) = each %aligners) {
        my $alignment_result_class = 'Genome::InstrumentData::AlignmentResult::' . Genome::InstrumentData::AlignmentResult->_resolve_subclass_name_for_aligner_name($aligner);
        next unless $alignment_result_class and ($alignment_result_class->can('prepare_reference_sequence_index') or $alignment_result_class->can('prepare_annotation_index'));

        while(my ($version, $references) = each %$versions) {
            while(my ($reference, $annotations) = each %$references) {
                while(my ($annotation, $params_lists) = each %$annotations) {
                    if($alignment_result_class->aligner_params_required_for_index){
                        for my $params (keys %$params_lists) {
                            my $step_data = $class->_generate_step($aligner, $version, $reference, $annotation, $params, $inputs);
                            $workflow_operations->{$step_data->{operation}} = $step_data;
                        }
                    } else {
                        my $step_data = $class->_generate_step($aligner, $version, $reference, $annotation, undef, $inputs);
                        $workflow_operations->{$step_data->{operation}} = $step_data;
                    }
                }
            }
        }
    }

    my $block_operation = $class->_wire_operations_to_master_workflow($master_workflow, $workflow_operations);

    return ($workflow_operations, $block_operation);
}

sub _generate_step {
    my $class = shift;
    my ($aligner, $version, $reference, $annotation, $params, $inputs) = @_;

    my $index_num = $class->next_counter_value;
    my $operation = Genome::WorkflowBuilder::Command->create(
        name => "$aligner index #" . $index_num,
        command => 'Genome::Model::ReferenceSequence::Command::CreateAlignerIndex',
    );

    #overwrite lsf_resource for star aligner
    if ($aligner eq 'star') {
        $operation->lsf_resource("-R \'select[ncpus>=12 && mem>=48000] span[hosts=1] rusage[mem=48000]\' -M 48000000 -n 12");
    }

    $operation->declare_constant(
        aligner_name => $aligner,
        aligner_version => $version,
        reference_sequence_build_id => $inputs->{$reference}->id,
        annotation_build_id => ($annotation? $inputs->{$annotation}->id : undef),
        aligner_params => (defined $params? $params : ''),
    );

    return {
        operation => $operation,
        index => $index_num,
    };
}

sub _wire_operations_to_master_workflow {
    my $class = shift;
    my $master_workflow = shift;
    my $index_operations = shift;

    my @block_operation_inputs = map { 'index_' . $_->{index} . '_result' } values %$index_operations;
    my $block_operation = Genome::WorkflowBuilder::Block->create(
        name => 'Wait for aligner indicies to be built',
        properties => [@block_operation_inputs,'force_fragment'],
    );
    $master_workflow->add_operation($block_operation);
    $master_workflow->connect_input(
        input_property => "m_force_fragment",
        destination => $block_operation,
        destination_property => "force_fragment",
    );

    for my $index_operation (values %$index_operations){
        $class->_wire_index_operation_to_master_workflow($master_workflow, $block_operation, $index_operation);
    }

    return $block_operation;
}

sub _wire_index_operation_to_master_workflow {
    my $class = shift;
    my $master_workflow = shift;
    my $block_operation = shift;
    my $operation = shift;

    $master_workflow->add_operation($operation->{operation});

    #result users are the same for all steps in workflow
    $master_workflow->connect_input(
        input_property => 'm_result_users',
        destination => $operation->{operation},
        destination_property => 'result_users',
    );

    $master_workflow->create_link(
        source => $operation->{operation},
        source_property => "result",
        destination => $block_operation,
        destination_property => join("_", "index", $operation->{index}, "result"),
    );

    return 1;
}

1;
