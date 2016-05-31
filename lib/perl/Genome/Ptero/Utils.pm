package Genome::Ptero::Utils;

use strict;
use warnings FATAL => 'all';
use Params::Validate qw(validate_pos :types);
use Ptero::Proxy::Workflow;
use Try::Tiny qw(try);
use JSON qw(from_json);

use Exporter 'import';

our @EXPORT_OK = qw(
    ptero_proxy
    ptero_workflow_url
    test_data_directory
    get_all_executions_for_proxy
    get_test_inputs
    get_test_outputs
    get_test_xml_filename
);

sub ptero_workflow_url {
    my ($workflow_name) = validate_pos(@_, {type => SCALAR});

    return sprintf("%s?name=%s", Genome::Config::get('ptero_workflow_submit_url'),
        $workflow_name);
}

sub ptero_proxy {
    my ($workflow_name) = validate_pos(@_, {type => SCALAR});

    my $proxy = try {
        return Ptero::Proxy::Workflow->new(
            ptero_workflow_url($workflow_name));
    };
    return $proxy;
}

sub get_all_executions_for_proxy {
    my ($workflow) = validate_pos(@_, {isa => 'Ptero::Proxy::Workflow'});

    my @executions = @{$workflow->workflow_executions};
    my @child_executions;

    for my $e (@executions) {
        for my $child (@{ $e->child_workflow_proxies }) {
            push @child_executions, get_all_executions_for_proxy($child);
        }
    }

    return (@executions, @child_executions);
}

sub get_test_xml_filename {
    my $name = shift;

    my $file = File::Spec->join(test_data_directory($name), 'workflow.xml');
    die "Cannot locate workflow.xml for workflow_test: $name" unless -e $file;
    return $file;
}

sub get_test_inputs {
    my $name = shift;
    my $file = File::Spec->join(test_data_directory($name), 'inputs.json');
    die "Cannot locate test inputs for workflow_test: $name" unless -e $file;
    return from_json(Genome::Sys->read_file($file));
}

sub get_test_outputs {
    my $name = shift;
    my $file = File::Spec->join(test_data_directory($name), 'outputs.json');
    die "Cannot locate test outputs for workflow_test: $name" unless -e $file;
    return from_json(Genome::Sys->read_file($file));
}

sub test_data_directory {
    my $name = shift;
    my $genome_dir = Genome->base_dir();

    return File::Spec->join($genome_dir, 'Ptero', 'workflow_tests', $name);
}


1;
