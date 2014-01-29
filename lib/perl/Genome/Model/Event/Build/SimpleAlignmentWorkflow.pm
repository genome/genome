package Genome::Model::Event::Build::SimpleAlignmentWorkflow;

use strict;
use warnings;
use Genome;

class Genome::Model::Event::Build::SimpleAlignmentWorkflow {
    is => ['Genome::Model::Event'],
};

sub help_brief {
    "Runs a simple alignment pipeline."
}

sub help_synopsis {
    return <<"EOS"
genome model build mymodel
EOS
}

sub help_detail {
    return <<"EOS"
One build of a simple alignment. 
EOS
}

sub execute {
    my $self = shift;

    $self->dump_status_messages(1);
    $self->dump_error_messages(1);

    # Verify the model
    my $model = $self->model;
    unless ($model) {
        $self->error_message("Failed to get a model for this build!");
        return;
    }
    my $build = $self->build;
    unless ($build) {
        $self->error_message("Failed to get a build object!");
        return;
    }

    my $data_directory = $self->build->data_directory;
    unless ($data_directory) {
        $self->error_message("Failed to get a data_directory for this build!");
        return;
    }

    # Get the processing profile and the params we care about
    my $processing_profile = $model->processing_profile;
    unless ($processing_profile) {
        $self->error_message("Failed to get a processing_profile object!");
        return;
    }

    my $workflow_log_directory = $self->build->workflow_log_directory_path;
    Genome::Sys->create_directory($workflow_log_directory);
    $self->debug_message("Created workflow log directory.");
    
    my $merged_file = $self->build->merged_alignment_file;

    my @instrument_data = $model->instrument_data;
    my @id_list;
    for my $instrument_data (@instrument_data) {
        $self->debug_message("Found instrument data id: ".$instrument_data->id);
        push(@id_list,$instrument_data->id); 
    }

    my $id_string = join(',',@id_list);
   
    $self->debug_message("Reference name: ".$model->reference_sequence_name );
    $self->debug_message("Inst. Id list: ".$id_string );
    $self->debug_message("Merged file: ".$merged_file );
    $self->debug_message("Workflow dir: ".$self->build->data_directory );
 
    my $workflow = Genome::Model::Tools::SimpleAlignment::PipelineInstrumentData->create(
        reference_name => $model->reference_sequence_name,
        instrument_data_id => $id_string,
        merged_file => $merged_file,
        working_directory =>  $self->build->data_directory 
    );

    unless ($workflow) {
        $self->error_message("Failed to create workflow!");
        return;
    }

    $workflow->execute();

    $self->debug_message("***Done with Build***");
    return 1;
}

1;
