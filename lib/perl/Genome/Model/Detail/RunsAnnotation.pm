package Genome::Model::Detail::RunsAnnotation;

use strict;
use warnings FATAL => 'all';
use Genome;
use Genome::Annotation::Dag qw(generate_dag);
use File::Basename qw(dirname);
use File::Slurp qw(write_file);

class Genome::Model::Detail::RunsAnnotation {
    is => 'Genome::ModelDeprecated',
    is_abstract => 1,

    has_param => [
        snvs_annotation_plan_name => {
            is => "Text",
            is_optional =>1,
            doc => "The name of the Annotation plan used to annotate and create reports for SNVs",
        },
        indels_annotation_plan_name => {
            is => "Text",
            is_optional =>1,
            doc => "The name of the Annotation plan used to annotate and create reports for INDELs",
        },
    ],
};

sub _get_annotation_plan {
    my ($self, $name) = @_;

    my $plan_file = File::Spec->join($self->_annotation_plan_search_dir, $name . '.yaml');
    unless (-f $plan_file) {
        die $self->error_message("Could not find annotation plan named ($name) in search directory (%s)", $self->_annotation_plan_search_dir);
    }
    my $plan = Genome::Annotation::Plan->create_from_file($plan_file);

    return $plan;
}

sub _annotation_plan_search_dir {
    my $self = shift;
    my $genome_base_dir = dirname(dirname(dirname(__FILE__)));
    return File::Spec->join($genome_base_dir, 'Annotation', 'plan_files');
}

sub annotation_plan {
    my ($self, $variant_type) = @_;
    return $self->_get_annotation_plan($self->annotation_plan_name($variant_type));
}

sub validated_annotation_plan {
    my ($self, $variant_type) = @_;
    my $plan = $self->annotation_plan($variant_type);
    $plan->validate();
    return $plan;
}

sub annotation_plan_name {
    my ($self, $variant_type) = @_;
    my $name_accessor = sprintf('%s_annotation_plan_name', $variant_type);
    return $self->$name_accessor;
}

sub workflow_xml_file {
    my $self = shift;
    my $basic_xml_file = shift;

    my $dag = Genome::WorkflowBuilder::DAG->from_xml_filename($basic_xml_file);

    $self->_connect_annotation_workflow($dag, 'snvs');
    $self->_connect_annotation_workflow($dag, 'indels');

    my $xml_file = Genome::Sys->create_temp_file_path;
    write_file($xml_file, $dag->get_xml);
    return $xml_file;
}

sub _connect_annotation_workflow {
    my ($self, $dag, $variant_type) = @_;

    if ($self->annotation_plan_name($variant_type)) {
        my $annotation_dag = generate_dag(
            $self->validated_annotation_plan($variant_type),
            $variant_type);

        $dag->add_operation($annotation_dag);

        my $dv_op = $dag->operation_named('Detect Variants');
        $dag->create_link(
            source => $dv_op,
            source_property => 'build_id',
            destination => $annotation_dag,
            destination_property => 'build_id',
        );

        for my $name qw(variant_type output_directory plan_json) {
            $dag->connect_input(
                input_property => sprintf('%s_%s', $variant_type, $name),
                destination => $annotation_dag,
                destination_property => $name,
            );
        }

        $dag->connect_output(
            output_property => sprintf('%s_output_directory', $variant_type),
            source => $annotation_dag,
            source_property => 'output_directory',
        );
    }
}

sub annotation_related_workflow_inputs {
    my ($self, $build) = @_;

    my %inputs;
    for my $variant_type qw(snvs indels) {
        if ($self->annotation_plan_name($variant_type)) {
            my $plan = $self->annotation_plan($variant_type);

            $inputs{sprintf('%s_variant_type', $variant_type)} = $variant_type;
            $inputs{sprintf('%s_output_directory', $variant_type)} =
                File::Spec->join($build->data_directory, 'reports', $variant_type);
            $inputs{sprintf('%s_plan_json', $variant_type)} = $plan->as_json();
        }
    }
    return %inputs;
}

sub __profile_errors__ {
    my $self = shift;
    my ($pp) = @_;

    my @errors = $self->SUPER::__profile_errors__(@_);
    for my $annotation_plan_name_accessor qw(snvs_annotation_plan_name indels_annotation_plan_name) {
        if (my $plan_name = $pp->$annotation_plan_name_accessor) {
            my $annotation_plan;
            eval {
                $annotation_plan = $self->_get_annotation_plan($plan_name);
            };
            if (my $error = $@) {
                push @errors, UR::Object::Tag->create(
                    type => 'error',
                    properties => [],
                    desc => $error,
                );
            } else {
                push @errors, $annotation_plan->__errors__;
            }
        }
    }

    return @errors;
}

1;
