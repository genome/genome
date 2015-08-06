package Genome::InstrumentData::Composite::Workflow::Generator::Base;

use strict;
use warnings;

use Genome;

class Genome::InstrumentData::Composite::Workflow::Generator::Base {
    is => 'UR::Singleton',
    is_abstract => 1,
};


my @_properties_for_add_link = ('id','destination','destination_property','source','source_property');
my @_properties_in_template_order;
my $_add_link_to_workflow_tmpl;
my $_add_link_workflow_id_pos;
my $_add_link_id_pos;
sub _add_link_to_workflow {
    my($class, $workflow_obj, %params) = @_;

    unless ($_add_link_to_workflow_tmpl) {
        my $tmpl = UR::BoolExpr::Template->resolve('Genome::WorkflowBuilder::Link', @_properties_for_add_link);
        unless ($tmpl) {
            Carp::croak('Cannot resolve BoolExpr template for adding Workflow::Link objects');
        }
        $tmpl = $tmpl->get_normalized_template_equivalent();
        for my $name ( @_properties_for_add_link ) {
            my $pos = $tmpl->value_position_for_property_name($name);
            if ($name eq 'id') {
                $_add_link_id_pos = $pos;
            }
            $_properties_in_template_order[$pos] = $name;
        }
        $_add_link_to_workflow_tmpl = $tmpl;
    }
    my @values = @params{@_properties_in_template_order};
    $values[$_add_link_id_pos] = UR::Object::Type->autogenerate_new_object_id_urinternal();
    my $rule = $_add_link_to_workflow_tmpl->get_rule_for_values(@values);
    return $workflow_obj->add_link(UR::Context->create_entity('Genome::WorkflowBuilder::Link', $rule));
}

my $counter = 0;
sub next_counter_value {
    return ++$counter;
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
