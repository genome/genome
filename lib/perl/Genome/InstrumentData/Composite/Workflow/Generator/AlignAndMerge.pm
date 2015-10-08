package Genome::InstrumentData::Composite::Workflow::Generator::AlignAndMerge;

use strict;
use warnings;
use Genome;
use POSIX qw(ceil);
use List::Util qw(sum);

class Genome::InstrumentData::Composite::Workflow::Generator::AlignAndMerge {
    is => 'Genome::InstrumentData::Composite::Workflow::Generator::Base',
};

sub generate {
    my $class = shift;
    my $tree = shift;
    my $input_data = shift;
    my $alignment_objects = shift;

    my $aligner_name = ucfirst($tree->{action}->[0]->{name});

    #Make a workflow with its input and output connectors
    my $input_properties = ['instrument_data', 'reference_sequence_build', $class->_general_workflow_input_properties];
    my $tree_properties = ['name', 'params', 'version'];
    my $workflow = Workflow::Model->create(
        name => $aligner_name,
        input_properties => [@$input_properties, @$tree_properties],
        optional_input_properties => $input_properties,
        output_properties => ['result_id'],
    );
    my $workflows = {};
    map { $workflows->{$_} = $workflow } @$alignment_objects;

    #Make a align_and_merge operation
    my $operation = $workflow->add_operation(
        name => "$aligner_name operation",
        operation_type => Workflow::OperationType::Command->create(
            command_class_name => 'Genome::InstrumentData::Command::AlignAndMerge',
        ),
    );

    my $instrument_data = $input_data->{instrument_data};
    my $lsf_resource_string = $class->_get_lsf_resource_string_for_aligner_and_instrument_data(
        $aligner_name,
        @$instrument_data
    );

    $operation->operation_type->lsf_resource($lsf_resource_string);

    #Connect input connectors to the operation
    my $inputs = [];
    for my $input_property (@$input_properties) {
        $workflow->add_link(
            left_operation => $workflow->get_input_connector,
            left_property => $input_property,
            right_operation => $operation,
            right_property => $input_property,
        );
        push @$inputs, ( 'm_' . $input_property => $input_data->{$input_property} );
    }
    for my $input_property (@$tree_properties) {
        $workflow->add_link(
            left_operation => $workflow->get_input_connector,
            left_property => $input_property,
            right_operation => $operation,
            right_property => $input_property,
        );
        push @$inputs, ( 'm_' . $input_property => $tree->{'action'}->[0]->{$input_property} );
    }

    #Connect output connectors to the operation
    $workflow->add_link(
        left_operation => $operation,
        left_property => 'result_id',
        right_operation => $workflow->get_output_connector,
        right_property => 'result_id',
    );

    if (exists $tree->{'action'}->[0]->{decoration}) {
        push @$inputs, Genome::InstrumentData::Composite::Decorator->decorate($operation, $tree->{'action'}->[0]->{decoration});
    }

    return $workflows, $inputs;
}

sub _get_lsf_resource_string_for_aligner_and_instrument_data {
    my $class = shift;
    my $aligner_name = shift;
    my @instrument_data = @_;

    my $merged_result_class = Genome::InstrumentData::Command::AlignAndMerge->merged_result_class($aligner_name);
    my $estimated_gtmp_bytes = sum(map { $merged_result_class->estimated_gtmp_for_instrument_data($_) } @instrument_data);
    return $class->_format_lsf_resource_string($estimated_gtmp_bytes);
}

sub _format_lsf_resource_string {
    my $class = shift;
    my $gtmp_bytes = shift;

    my $cpus = 8;
    my $mem_gb = 60;
    my $queue = Genome::Config::get('lsf_queue_alignment_default');

    my $gtmp_kb = ceil($gtmp_bytes / 1024);
    my $gtmp_mb = ceil($gtmp_kb / 1024);
    my $gtmp_gb = ceil($gtmp_mb / 1024);

    my $mem_mb = $mem_gb * 1024;
    my $mem_kb = $mem_mb * 1024;

    my $select  = "select[ncpus >= $cpus && mem >= $mem_mb && gtmp >= $gtmp_gb] span[hosts=1]";
    my $rusage  = "rusage[mem=$mem_mb, gtmp=$gtmp_gb]";
    my $options = "-M $mem_kb -n $cpus -q $queue";

    my $required_usage = "-R \'$select $rusage\' $options";

    return $required_usage;
}

# sub _workflow_

1;
