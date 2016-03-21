package Genome::InstrumentData::Composite::Workflow::Generator::Merge;

use strict;
use warnings;

use Genome;

class Genome::InstrumentData::Composite::Workflow::Generator::Merge {
    is => 'Genome::InstrumentData::Composite::Workflow::Generator::Base',
};

sub generate {
    my $class = shift;
    my $master_workflow = shift;
    my $tree = shift;
    my $group = shift;
    my $alignment_objects = shift;

    my $merge_operation;
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
            }
        }

        $merge_operation = $class->_generate_merge_operation($merge_tree, $group);
        $class->_wire_merge_operation_to_master_workflow($master_workflow, $merge_operation);
    }

    return ($merge_operation, \@inputs);
}

my $_generate_merge_operation_tmpl;
sub _generate_merge_operation {
    my $class = shift;
    my $merge_tree = shift;
    my $grouping = shift;

    unless ($_generate_merge_operation_tmpl) {
        $_generate_merge_operation_tmpl = UR::BoolExpr::Template->resolve('Genome::WorkflowBuilder::Command', 'id','name','command')->get_normalized_template_equivalent();
    }

    my $operation = Genome::WorkflowBuilder::Command->create(
        $_generate_merge_operation_tmpl->get_rule_for_values(
                'Genome::InstrumentData::Command::MergeAlignments',
                UR::Object::Type->autogenerate_new_object_id_urinternal(),
                'merge_' . $grouping,
        )
    );

    return $operation;
}

sub _wire_merge_operation_to_master_workflow {
    my $class = shift;
    my $master_workflow = shift;
    my $merge = shift;

    $master_workflow->add_operation($merge);
    for my $property ($class->_merge_workflow_input_properties) {
        $master_workflow->connect_input(
            input_property => 'm_' . $property,
            destination => $merge,
            destination_property => $property,
            is_optional => $class->_is_duplication_parameter($property),
        );
    }

    return 1;
}

sub _is_duplication_parameter {
    my $class = shift;
    my $param = shift;

    return (index($param, 'duplication_handler') == 0);
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

1;
