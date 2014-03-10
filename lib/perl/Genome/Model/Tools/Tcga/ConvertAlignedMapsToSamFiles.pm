package Genome::Model::Tools::Tcga::ConvertAlignedMapsToSamFiles;

use strict;
use warnings;

use Genome;
use Workflow;
use IO::File;

class Genome::Model::Tools::Tcga::ConvertAlignedMapsToSamFiles {
    is  => ['Command'],
    has => [
        model_id => {
            is  => 'String',
            is_input => '1',
            doc => 'The model id.',
        },
        working_directory => {
            is => 'String',
            is_input => '1',
            doc => 'The working directory.',
        },
        aligned_sam_file_directory => {
            is => 'String',
            is_output =>1,
            is_optional =>1,
            doc => 'The directory where all the resulting sam files will be generated.', 
        }, 
    ],
};

sub help_brief {
    'Convert Maq map files into the TCGA format.';
}

sub help_detail {
    return <<EOS
    Convert Maq map files into the TCGA format.
EOS
}


sub execute {
    my $self = shift;

    $self->dump_status_messages(1);
    my $model_id = $self->model_id;
    my $model = Genome::Model->get($model_id);
    die "Model $model_id is not defined. Quitting." unless defined($model);
     
    my @instrument_data = $model->instrument_data;
    $self->debug_message("There are " . scalar(@instrument_data) . " id assignemnts for model id $model_id\n");
    my $build = $model->last_complete_build;
    unless ($build) {
        die "No successful build of model " . $model->__display_name;
    }
    my $count=0;

    my @alignments;
    for my $instrument_data (@instrument_data) {
        my $idid = $instrument_data->id;
        my $alignment = $build->alignment_results_for_instrument_data($instrument_data);
        my $alignment_directory = $alignment->output_dir;
        
        push (@alignments, "$idid|$alignment_directory");
    }

    $self->debug_message("Alignment info sent to workers: ".join("\n",@alignments));
    $self->debug_message("Working dir sent to workers: ".$self->working_directory);

    require Workflow::Simple;

    my $op = Workflow::Operation->create(
        name => 'Generate per lane sams',
        operation_type => Workflow::OperationType::Command->get('Genome::Model::Tools::Tcga::ConvertAlignedMapsToSamFilesWorker')
    );

    $op->parallel_by('alignment_info');

    my $output = Workflow::Simple::run_workflow_lsf(
        $op,
        'alignment_info'  => \@alignments,
        'working_directory' => $self->working_directory, 
    );

    #check workflow for errors 
    if (!defined $output) {
        foreach my $error (@Workflow::Simple::ERROR) {
            $self->error_message($error->error);
        }
        return;
    } else {
        $self->debug_message("Workflow completed with no errors.");
    }


    $self->aligned_sam_file_directory($self->working_directory."/aligned/");

    return 1;

}
1;
