package Genome::InstrumentData::Composite::Workflow::Generator::Base;

use strict;
use warnings;

use Genome;

class Genome::InstrumentData::Composite::Workflow::Generator::Base {
    is => 'UR::Singleton',
    is_abstract => 1,
};


my @_properties_for_add_link = ('id','workflow_model_id','left_workflow_operation_id','left_property',
                                'right_workflow_operation_id','right_property');
my @_properties_in_template_order;
my $_add_link_to_workflow_tmpl;
my $_add_link_workflow_id_pos;
my $_add_link_id_pos;
sub _add_link_to_workflow {
    my($class, $workflow_obj, %params) = @_;

    unless ($_add_link_to_workflow_tmpl) {
        my $tmpl = UR::BoolExpr::Template->resolve('Workflow::Link', @_properties_for_add_link);
        unless ($tmpl) {
            Carp::croak('Cannot resolve BoolExpr template for adding Workflow::Link objects');
        }
        $tmpl = $tmpl->get_normalized_template_equivalent();
        for my $name ( @_properties_for_add_link ) {
            my $pos = $tmpl->value_position_for_property_name($name);
            if ($name eq 'workflow_model_id' ) {
                $_add_link_workflow_id_pos = $pos;
            } elsif ($name eq 'id') {
                $_add_link_id_pos = $pos;
            }
            $_properties_in_template_order[$pos] = $name;
        }
        $_add_link_to_workflow_tmpl = $tmpl;
    }
    my @values = @params{@_properties_in_template_order};
    $values[$_add_link_workflow_id_pos] = $workflow_obj->id;
    $values[$_add_link_id_pos] = UR::Object::Type->autogenerate_new_object_id_urinternal();
    my $rule = $_add_link_to_workflow_tmpl->get_rule_for_values(@values);
    return UR::Context->create_entity('Workflow::Link', $rule);
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
    if($@){
        $class->error_message("Could not create an instance of " . $class_name);
        die $class->error_message;
    }

    unless ($_generate_operation_tmpl) {
        # 'id' will still be first, since they're sorted alpha order
#ccc
        $_generate_operation_tmpl
            = UR::BoolExpr::Template->resolve('Workflow::Operation', 'id','name','workflow_operationtype_id')->get_normalized_template_equivalent();
    }

    my $operation = Workflow::Operation->create(
        $_generate_operation_tmpl->get_rule_for_values(UR::Object::Type->autogenerate_new_object_id_urinternal(),
        $class->get_unique_action_name($action),
        Workflow::OperationType::Command->get($class_name)->id),
    );

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

my $counter = 0;
sub next_counter_value {
    return ++$counter;
}

sub get_unique_action_name {
    my $class = shift;
    my $action = shift;

    my $index = $class->next_counter_value;

    my $name = $action->{name};
    return join('_', $name, $index);
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

sub _general_workflow_input_properties {
    my $class = shift;

    return qw(
        picard_version
        samtools_version
        bedtools_version
        force_fragment
        trimmer_name
        trimmer_version
        trimmer_params
        result_users
    );
}

my $REFINEMENT_INPUT_PROPERTY_SEPARATOR = ':';
sub _construct_refiner_input_property {
    my $class = shift;
    my $property = shift;
    my $refiner = shift;

    return join($REFINEMENT_INPUT_PROPERTY_SEPARATOR, $property, $refiner);
}

1;
