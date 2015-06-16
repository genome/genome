package Genome::WorkflowBuilder::Block;

use strict;
use warnings;

use Genome;

class Genome::WorkflowBuilder::Block {
    is => 'Genome::WorkflowBuilder::Converge',
};

sub get_ptero_builder_task {
    require Ptero::Builder::Detail::Workflow::Task;

    my $self = shift;

    $self->validate;

    my %params = (
        name => $self->name,
        methods => [
            $self->_get_ptero_block_method(),
        ],
    );
    if (defined $self->parallel_by) {
        $params{parallel_by} = $self->parallel_by;
    }
    return Ptero::Builder::Detail::Workflow::Task->new(%params);
}

sub _get_ptero_block_method {
    require Ptero::Builder::Detail::Workflow::Converge;

    my $self = shift;
    my ($output_name) = $self->output_properties;
    return Ptero::Builder::Detail::Workflow::Block->new(
        name => 'block',
        parameters => {
            input_names => [$self->input_properties],
            output_name => $output_name,
        },
    );
}

1;
