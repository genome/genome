package Genome::Model::Tools::Tcga::MergePerLaneSamFiles;

use strict;
use warnings;

use Genome;
use Workflow;
use IO::File;

class Genome::Model::Tools::Tcga::MergePerLaneSamFiles {
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
            is_input =>1,
            doc => 'The directory where all the resulting sam files will be generated.', 
        }, 
        unaligned_sam_file_directory => {
            is => 'String',
            is_input =>1,
            doc => 'The directory where all the resulting sam files will be generated.', 
        },
        read_group_directory => {
            is => 'String',
            is_input =>1,
            doc => 'All of the readgroup entries for the TCGA header.',
        },
        program_group_directory => {
            is => 'String',
            is_input =>1,
            doc => 'All of the programgroup entries for the TCGA header.',
        },
        final_file => {
            is => 'String',
            is_output =>1,
            is_optional =>1,
            doc => 'All of the required components merged together. ',
        },

    ],
    has_param => [
           lsf_resource => {
            default_value => 'select[model!=Opteron250 && type==LINUX64] rusage[mem=4000]',
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
    my $model = Genome::Model::ReferenceAlignment->get($model_id);
    die "Reference alignment model $model_id is not defined. Quitting." unless defined($model);
    
    my $seq_dict_sam_file = $model->reference_sequence_build->get_sequence_dictionary("sam");
 
    my @instrument_data = $model->instrument_data;
    $self->status_message("There are " . scalar(@instrument_data) . " id assignemnts for model id $model_id\n");

    my $build = $model->last_complete_build;
    unless ($build) {
        die "No sucessful build of model " . $model->__display_name__;
    }

    my @ids;
    for my $instrument_data (@instrument_data) {
        my $alignment = $build->alignment_results_for_instrument_data;
        next unless $alignment;
        push @ids, $instrument_data->id;
    }

    $self->status_message("Beginning per lane merge of the seq dict sam file, the rg rile, pg file, aligned file and unaligned file.");
    $self->status_message("Per lane seq id's to merge: " . join("\n",@ids));
    $self->status_message("Working dir sent to workers: " . $self->working_directory);

    require Workflow::Simple;

    my $op = Workflow::Operation->create(
        name => 'Generate per lane sams',
        operation_type => Workflow::OperationType::Command->get('Genome::Model::Tools::Tcga::MergePerLaneSamFilesWorker')
    );

    $op->parallel_by('seq_ids');

    my $output = Workflow::Simple::run_workflow_lsf(
        $op,
        'seq_ids'  => \@ids,
        'working_directory' => $self->working_directory, 
        'seq_dict_sam_file' => $seq_dict_sam_file, 
    );

    #check workflow for errors 
    if (!defined $output) {
        foreach my $error (@Workflow::Simple::ERROR) {
            $self->error_message($error->error);
        }
        return;
    } else {
        $self->status_message("Workflow completed with no errors.");
    }

    #merge all the bam files

    my $fixmated_file = $self->working_directory."/tcga_fixmated.bam"; 
    my $bam_file_dir = $self->working_directory."/bams/";
    if (!-s $fixmated_file) {
        my @bam_files = <$bam_file_dir/*.bam>;
        my $bam_file_list = join(" I=",@bam_files);
        my $merge_cmd = "java -Xmx2g -cp /gsc/scripts/lib/java/samtools/picard-tools-1.04/MergeSamFiles.jar net.sf.picard.sam.MergeSamFiles MSD=true SO=coordinate AS=true VALIDATION_STRINGENCY=SILENT O=$fixmated_file I=$bam_file_list";

        $self->status_message("Merging per lane bams into a single, large fixmated file with: $merge_cmd");
        my $merge_rv = Genome::Sys->shellcmd(cmd=>$merge_cmd,input_files=>\@bam_files,output_files=>[$fixmated_file]);
        if ($merge_rv != 1) {
            die "Merge with command: $merge_cmd failed with rv $merge_rv";
        } else {
            $self->status_message("Merge succeeded. $fixmated_file created.");
        } 

    } else {
        $self->status_message("The bam file: $fixmated_file already exists.  Skipping the generation of this file.");
    }

    my $markdup_file = $self->working_directory."/tcga_markdup.bam";   
    my $metrics_file = $self->working_directory."/markdup.metrics";
    my $log_file = $self->working_directory."/markdup.log";

    $self->status_message("Beginning MarkDuplicates.  Attempting to generate $markdup_file.");
    if (!-s $markdup_file) { 
        my $markdup_cmd = Genome::Model::Tools::Sam::MarkDuplicates->create(file_to_mark=>$fixmated_file,
            marked_file=>$markdup_file,
            metrics_file=>$metrics_file,
            remove_duplicates=>0,
            max_jvm_heap_size=>2,        
            tmp_dir=>$self->working_directory."/tmp/", 
        );

        my $markdup_rv = $markdup_cmd->execute();
        if ($markdup_rv ne 1) {
            $self->status_message("MarkDuplicates Error:  The return value was: ".$markdup_rv);
            die;
        } else {
            $self->status_message("MarkDuplicates succeeded! Return value:".$markdup_rv);
        }  

    } else {
        $self->status_message("The bam file: $markdup_file already exists.  Skipping the generation of this file.");
    }

    $self->final_file($markdup_file);
    return 1;

}
1;
