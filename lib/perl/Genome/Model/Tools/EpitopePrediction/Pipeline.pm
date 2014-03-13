package Genome::Model::Tools::EpitopePrediction::Pipeline;

use strict;
use warnings;

use Genome;
use Workflow::Simple;

class Genome::Model::Tools::EpitopePrediction::Pipeline {
    is => 'Command::V2',
    doc => 'Run the epitope binding prediction pipeline',
    has => [
        output_directory => {
            is => 'Text',
            is_optional => 0,
            doc => 'the directory where you want results stored',
        },
    ],
};

sub generate_result {
    my ($self) = @_;

    $self->debug_message("Validating Inputs...");
    $self->_validate_inputs();

    $self->debug_message("Constructing Workflow...");
    my $workflow = $self->_construct_workflow();

    $self->debug_message("Getting Workflow Inputs...");
    my $inputs = $self->_get_workflow_inputs();

    $self->debug_message("Running Workflow...");
    my $result = Workflow::Simple::run_workflow_lsf($workflow, %$inputs);

    unless($result){
        $self->error_message( join("\n", map($_->name . ': ' . $_->error, @Workflow::Simple::ERROR)) );
        die $self->error_message("Workflow did not return correctly.");
    }

    return 1;
}

sub _construct_workflow {
    my ($self) = @_;

    my $xml = __FILE__ . '.xml';
    my $workflow = Workflow::Operation->create_from_xml($xml);
    $workflow->log_dir($self->output_directory);

    return $workflow;
}

sub _validate_inputs {
    my $self = shift;
    die $self->error_message('TODO: Validate inputs if needed');
}

sub _get_workflow_inputs {
    my $self = shift;
    die $self->error_message('TODO: Fill in the input hash for the workflow -- all input connector properties must be set');
    my %inputs = (
    );

    return \%inputs;
}

1;
