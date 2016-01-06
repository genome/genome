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

    my $master_workflow = Genome::WorkflowBuilder::DAG->create(
        name => 'Master Alignment Dispatcher',
    );

    my $block_operation = Genome::InstrumentData::Composite::Workflow::Generator::AlignerIndex->generate($master_workflow, $tree, $input_data);

    my %inputs;

    my $objects_by_group = $class->_alignment_objects($tree, $input_data, $merge_group);
    for my $group (keys %$objects_by_group) {
        my $alignment_objects = $objects_by_group->{$group};

        my $alignment_generator = 'Genome::InstrumentData::Composite::Workflow::Generator::' . string_to_camel_case($tree->{action}->[0]->{type});
        my ($next_object_workflows, $next_object_inputs) = $alignment_generator->generate($master_workflow, $block_operation, $tree, $input_data, $alignment_objects);
        %inputs = (%inputs, @$next_object_inputs);

        my ($next_merge_operation, $next_merge_inputs) = Genome::InstrumentData::Composite::Workflow::Generator::Merge->generate($master_workflow, $tree, $group, $alignment_objects);
        if($next_merge_operation) {
            %inputs = (%inputs, @$next_merge_inputs);

            $class->_wire_object_workflows_to_merge_operations($master_workflow, $next_object_workflows, $next_merge_operation, $alignment_objects, $group);

            my ($next_refinement_operations, $next_refinement_inputs, $next_refiners) = Genome::InstrumentData::Composite::Workflow::Generator::Refine->generate($master_workflow, $tree, $input_data, $group, $alignment_objects);
            if(@$next_refinement_operations) {
                %inputs = (%inputs, @$next_refinement_inputs);

                $class->_wire_merge_to_refinement_operation($master_workflow, $next_merge_operation, $next_refinement_operations->[0]);
            } else {
                $class->_wire_merge_to_output($master_workflow, $next_merge_operation);
            }
        }
    }

    return ($master_workflow, $class->_resolve_inputs(\%inputs, $tree->{api_version}, $input_data));
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

sub _resolve_inputs {
    my $class = shift;
    my $inputs = shift;
    my $api_version = shift;
    my $input_data = shift;

    #add the global inputs
    my $api_inputs = $class->inputs_for_api_version($api_version);
    my %inputs = (%$api_inputs, %$input_data);
    for my $input ($class->_general_workflow_input_properties) {
            $inputs->{'m_' . $input} = $inputs{$input};
    }

    return $inputs;
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

sub _wire_merge_to_refinement_operation {
    my $class = shift;
    my $master_workflow = shift;
    my $merge_operation = shift;
    my $refinement_operation = shift;

    $master_workflow->create_link(
        source => $merge_operation,
        source_property => 'result_id',
        destination => $refinement_operation,
        destination_property => 'input_result_id',
    );

    return 1;
}

sub _wire_merge_to_output {
    my $class = shift;
    my $master_workflow = shift;
    my $merge_operation = shift;

    for my $property ($merge_operation->output_properties) {
        $master_workflow->connect_output(
            source => $merge_operation,
            source_property => $property,
            output_property => 'm_' . join('_', $property, $merge_operation->name),
        );
    }
}

sub _wire_object_workflows_to_merge_operations {
    my $class = shift;
    my $master_workflow = shift;
    my $object_workflows = shift;
    my $merge_operation = shift;
    my $alignment_objects = shift;
    my $group = shift;

    my @converge_inputs;
    for my $o (@$alignment_objects) {
        my $align_wf = $object_workflows->{$o};

        push @converge_inputs, $align_wf->output_properties;
    }

    my $converge_operation = Genome::WorkflowBuilder::Converge->create(
        name => join('_', 'converge', $group),
        input_properties => \@converge_inputs,
        output_properties => ['alignment_result_ids'],
    );
    $master_workflow->add_operation($converge_operation);
    $master_workflow->create_link(
        source => $converge_operation,
        source_property => 'alignment_result_ids',
        destination => $merge_operation,
        destination_property => 'alignment_result_ids',
    );

    for my $o (@$alignment_objects) {
        my $align_wf = $object_workflows->{$o};

        for my $property ($align_wf->output_properties) {
            $master_workflow->create_link(
                source => $align_wf,
                source_property => $property,
                destination => $converge_operation,
                destination_property => $property,
            );
        }
    }

    return 1;
}

1;
