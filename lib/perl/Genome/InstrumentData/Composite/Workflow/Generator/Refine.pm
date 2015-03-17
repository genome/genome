package Genome::InstrumentData::Composite::Workflow::Generator::Refine;

use strict;
use warnings;

use Genome;

class Genome::InstrumentData::Composite::Workflow::Generator::Refine {
    is => 'Genome::InstrumentData::Composite::Workflow::Generator::Base',
};

sub generate {
    my $class = shift;
    my $tree = shift;
    my $input_data = shift;
    my $group = shift;
    my $alignment_objects = shift;

    my @refinement_operations;
    my @inputs;

    my @refiners;

    if(exists $tree->{then}) {
        my $merge_tree = $tree->{then};
        my $next_op = $merge_tree;
        while(exists $next_op->{then}) {
            $next_op = $next_op->{then};

            if($next_op->{type} eq 'refine') {
                my $refiner = $next_op->{name};
                if (grep { $_ eq $refiner } @refiners) {
                    die $class->error_message("Refiner $refiner is used more than once");
                }
                push @refiners, $refiner;
                push @inputs, (
                    $class->_construct_refiner_input_property('m_refiner_name', $refiner)   => $next_op->{name},
                    $class->_construct_refiner_input_property('m_refiner_params', $refiner)  => $next_op->{params},
                    $class->_construct_refiner_input_property('m_refiner_version', $refiner) => $next_op->{version},
                    $class->_construct_refiner_input_property('m_refiner_known_sites_ids', $refiner) => [ map { $_->id } @{$input_data->{$next_op->{known_sites}}} ],
                );

                my $key = $class->refiner_key($group, $refiner);
                push @refinement_operations, $class->_generate_refinement_operation($merge_tree, $key);
            }
        }
    }

    return (\@refinement_operations, \@inputs, \@refiners);
}

my $_refinealignments_command_id;
my $_generate_refinement_operation_tmpl;
sub _generate_refinement_operation {
    my $class = shift;
    my $merge_tree = shift;
    my $refiner = shift;

    unless ($_generate_refinement_operation_tmpl) {
        $_generate_refinement_operation_tmpl = UR::BoolExpr::Template->resolve('Workflow::Operation', 'id','name','workflow_operationtype_id')->get_normalized_template_equivalent();
        $_refinealignments_command_id = Workflow::OperationType::Command->get('Genome::InstrumentData::Command::RefineReads')->id;
    }

    my $operation = Workflow::Operation->create(
        $_generate_refinement_operation_tmpl->get_rule_for_values(
                UR::Object::Type->autogenerate_new_object_id_urinternal(),
                'refine_' . $refiner,
                $_refinealignments_command_id
        )
    );

    return $operation;
}

1;
