package Genome::InstrumentData::Composite::Workflow::Generator::Merge;

use strict;
use warnings;

use Genome;

class Genome::InstrumentData::Composite::Workflow::Generator::Merge {
    is => 'Genome::InstrumentData::Composite::Workflow::Generator::Base',
};

sub generate {
    my $class = shift;
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
    }

    return ($merge_operation, \@inputs);
}

my $_mergealignments_command_id;
my $_generate_merge_operation_tmpl;
sub _generate_merge_operation {
    my $class = shift;
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

1;
