package Genome::InstrumentData::Composite::Workflow::Generator::Refine;

use strict;
use warnings;

use Genome;

class Genome::InstrumentData::Composite::Workflow::Generator::Refine {
    is => 'Genome::InstrumentData::Composite::Workflow::Generator::Base',
};

sub generate {
    my $class = shift;
    my $master_workflow = shift;
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

    if(@refinement_operations) {
        $class->_wire_refinement_operation_to_master_workflow($master_workflow, \@refinement_operations, \@refiners);
        $class->_wire_refinement_to_refinement_operations($master_workflow, \@refinement_operations);
    }

    return (\@refinement_operations, \@inputs, \@refiners);
}

my $_generate_refinement_operation_tmpl;
sub _generate_refinement_operation {
    my $class = shift;
    my $merge_tree = shift;
    my $refiner = shift;

    unless ($_generate_refinement_operation_tmpl) {
        $_generate_refinement_operation_tmpl = UR::BoolExpr::Template->resolve('Genome::WorkflowBuilder::Command', 'id','name','command')->get_normalized_template_equivalent();
    }

    my $operation = Genome::WorkflowBuilder::Command->create(
        $_generate_refinement_operation_tmpl->get_rule_for_values(
                'Genome::InstrumentData::Command::RefineReads',
                UR::Object::Type->autogenerate_new_object_id_urinternal(),
                'refine_' . $refiner,
        )
    );

    return $operation;
}

sub refiner_key {
    my ($class, $group, $refiner_name) = @_;
    return join ("_", $group, $refiner_name);
}

sub _wire_refinement_operation_to_master_workflow {
    my $class = shift;
    my $master_workflow = shift;
    my $refinement_operation = shift;
    my $refiners = shift;

    for my $refinement (@$refinement_operation) {
        $master_workflow->add_operation($refinement);
        my ($refiner) = grep { $refinement->name =~ /$_/ } @$refiners;

        for my $property ($class->_base_refinement_workflow_input_properties) {
            my $source_property = "m_" . $class->_construct_refiner_input_property($property, $refiner);
            my $destination_property = $property;
            $class->_add_link_to_workflow($master_workflow,
                source_property => $source_property,
                destination => $refinement,
                destination_property => $destination_property,
            );
        }

        $class->_add_link_to_workflow($master_workflow,
            source_property => 'm_result_users',
            destination => $refinement,
            destination_property => 'result_users',
        );
    }

    my $last_refinement = $refinement_operation->[-1];
    for my $property ($last_refinement->output_properties) {
        $class->_add_link_to_workflow($master_workflow,
            source => $last_refinement,
            source_property => $property,
            destination_property => 'm_' . join('_', $property, $last_refinement->name),
        );
    }

    return 1;
}

sub _wire_refinement_to_refinement_operations {
    my $class = shift;
    my $master_workflow = shift;
    my $refinement_operations = shift;

    if (scalar @$refinement_operations > 1) {
        for my $i (1..$#$refinement_operations) {
            $class->_add_link_to_workflow($master_workflow,
                source => $refinement_operations->[$i-1],
                source_property => 'result_id',
                destination => $refinement_operations->[$i],
                destination_property => 'input_result_id',
            );
        }
    }

    return 1;
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

my $REFINEMENT_INPUT_PROPERTY_SEPARATOR = ':';
sub _construct_refiner_input_property {
    my $class = shift;
    my $property = shift;
    my $refiner = shift;

    return join($REFINEMENT_INPUT_PROPERTY_SEPARATOR, $property, $refiner);
}

1;
