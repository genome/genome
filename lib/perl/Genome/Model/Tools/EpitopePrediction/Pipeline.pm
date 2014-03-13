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
            doc => 'the directory where you want results stored',
        },
    ],
};

sub execute {
    my $self = shift;

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

    unless (-s $self->input_tsv_file) {
        die $self->error_message("Input tsv file %s does not exist or has no size", $self->input_tsv_file);
    }

    unless (Genome::Sys->create_directory($self->output_directory)) {
        die $self->error_message("Coult not create directory (%s)", $self->output_directory);
    }

    # TODO make sure anno db makes sense
    # TODO make sure anno db version makes sense
    # TODO make sure length makes sense
    # TODO make sure allele makes sense
    # TODO make sure epitope_length makes sense
    # TODO make sure netmhc_version makes sense
    die $self->error_message('TODO: Validate inputs');

    return 1;
}

sub _get_workflow_inputs {
    my $self = shift;

    my %inputs = (
        input_tsv_file => $self->input_tsv_file,
        output_directory => $self->output_directory,
        anno_db => $self->anno_db,
        anno_db_version => $self->anno_db_version,
        length => $self->length,
        allele => $self->allele,
        epitope_length => $self->epitope_length,
        netmhc_version => $self->netmhc_version,
    );

    return \%inputs;
}

1;
