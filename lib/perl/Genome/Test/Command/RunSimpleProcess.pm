package Genome::Test::Command::RunSimpleProcess;

use strict;
use warnings FATAL => 'all';
use Genome;
use Genome::WorkflowBuilder::DAG;
use Genome::WorkflowBuilder::Command;

class Genome::Test::Command::RunSimpleProcess {
    is => 'Command::V2',

    has_input => [
        statement => {
            is => 'Text',
            doc => "Something to expect in the workflow outputs",
        },
        wait_and_check_outputs => {
            is => 'Boolean',
            default => '0',
        },
        fail_intentionally => {
            is => 'Boolean',
            default => '0',
        }
    ],
    doc => 'Runs a simple Genome::Process',
};

sub execute {
    my $self = shift;

    my $p = Genome::Process::Test::Process->create(
        statement => $self->statement,
    );

    $self->status_message("Constructing workflow from inputs. (this may take a while...)");
    my %run_params = (
        workflow_xml => $self->dag->get_xml,
        workflow_inputs => {statement => $self->statement,
                            fail_intentionally => $self->fail_intentionally},
    );
    if ($self->wait_and_check_outputs) {
        my $outputs = $p->run_and_wait(%run_params);
        die "Unexpected outputs!" unless $outputs->{statement} eq $self->statement;
    } else {
        $p->run(%run_params);
    }

    return $p;
}

sub dag {
    my $self = shift;

    my $dag = Genome::WorkflowBuilder::DAG->create(
        name => 'SimpleProcess',
    );

    my $command = Genome::WorkflowBuilder::Command->create(
        name => 'Run Echo Command',
        command => 'Genome::Test::Command::Echo',
    );
    $dag->add_operation($command);

    $dag->connect_input(
        input_property => 'statement',
        destination => $command,
        destination_property => 'input_statement',
    );
    $dag->connect_input(
        input_property => 'fail_intentionally',
        destination => $command,
        destination_property => 'fail_intentionally',
    );

    $dag->connect_output(
        output_property => 'statement',
        source => $command,
        source_property => 'output_statement',
    );

    return $dag;
}

1;
