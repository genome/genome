package Genome::Model::Tools::Cmds::CompileCnaOutput;

use warnings;
use strict;
use Genome;
use Cwd;

class Genome::Model::Tools::Cmds::CompileCnaOutput {
    is => 'Genome::Command::Base',
    has => [
    window_size=> {
        type => 'String',
        is_optional => 1,
        doc => "The window size is how far away one sample position is picked",
        default_value => 10000,
    },
    output_dir=> {
        type => 'String',
        is_input => 1,
        is_output => 1,
        is_optional => 1,
        default => getcwd(),
        doc => "Where the bam-to-cna file will be created, when window-size is not 10K, default is current working directory",
    },
    regenerate_missing_output => {
        type => 'Boolean',
        is_optional => 1,
        default => 0,
        doc => "If set to true, any missing copy number output files will be regenerated.",
    },
    force => {
        type => 'Boolean',
        is_optional => 1,
        is_input => 1,
        default => 0,
        doc => "If --force is set, continue despite missing somatic models or builds",
    },
    ],
    has_many => [
    models => { 
        shell_args_position => 1, 
        is => 'Genome::Model::Somatic', 
        is_input => 1 
    },
    ]
};

sub help_brief {
    'Create links to bam-to-cna files from somatic builds.'
}

sub help_detail {
    'This script checks for the copy number output and the associated .png image in somatic builds, and then generates a symbolic link to them if they exist. If they do not exist, the script will submit jobs to run the bam-to-cna script. Messages will indicate the status of the builds and actions taken, and the std_out of the submitted jobs will be printed in the working directory, so that you will be able to see their status and somatic model_ids associated. You will want to run this script again to create the missing symbolic links when all of your builds have the appropriate copy_number_output files.'
}

sub execute {
    my $self = shift;
    my @somatic_models;

    my $window_size = $self->window_size;
    my $output_dir = $self->output_dir;

    for my $somatic_model ($self->models) {
        my $somatic_model_id = $somatic_model->id;
        my $somatic_build = $somatic_model->last_succeeded_build;
        unless ($somatic_build) {
            $self->error_message("No succeeded build found for somatic model id $somatic_model_id.");
            if ($self->force) {
                $self->debug_message("Continuing despite the error above, skipping this model");
                next;
            } else {
                die;
            }
        }
        my $somatic_build_id = $somatic_build->id or die "No build id found in somatic build object for somatic model id $somatic_model_id.\n";
        $self->debug_message("Last succeeded build for somatic model $somatic_model_id is build $somatic_build_id. ");

        ## must changed in the new version
        my $cn_data_file = $somatic_build->somatic_workflow_input("copy_number_output") or die "Could not query somatic build for copy number output.\n";
        my $cn_png_file = $cn_data_file.".png";
        print "$cn_data_file\n";
        
        #if files are found (bam-to-cna has been run correctly already), create link to the data in current folder
        if (-s $cn_data_file && -s $cn_png_file) {
            if ($window_size == 10000){       
                my $link_name = "$output_dir/$somatic_model_id.copy_number.csv";
                if (-e $link_name) {
                    $self->debug_message("Link $link_name already found in dir.\n");
                } else {
                    `ln -s $cn_data_file $link_name`;
                    $self->debug_message("Link to copy_number_output created.\n");
                }
            } else{
                $self->debug_message("Window size is not 10000... this necessitates regeneration of copy number output for somatic model id $somatic_model_id");
                if($self->regenerate_missing_output) {
                    $self->regenerate_cnv_output($somatic_build);
                }
            }
        } else{
            $self->debug_message("Copy number output not found for build $somatic_build_id. ");
            if($self->regenerate_missing_output) {
                $self->regenerate_cnv_output($somatic_build);
            }
        }

    }   
    return 1;
}

# Generates the copy number output for a set of bams in the output dir (instead of symlinking existing output)
sub regenerate_cnv_output {
    my $self = shift;
    my $somatic_build_object = shift;
    my $output_dir = $self->output_dir;

    #get tumor and normal bam files
    $self->debug_message("Build new copy number output file in $output_dir. ");
    my $tumor_build = $somatic_build_object->tumor_build or die "Cannot find tumor model.\n";
    my $normal_build = $somatic_build_object->normal_build or die "Cannot find normal model.\n";
    my $tumor_bam = $tumor_build->whole_rmdup_bam_file or die "Cannot find tumor .bam.\n";
    my $normal_bam = $normal_build->whole_rmdup_bam_file or die "Cannot find normal .bam.\n";
    my $somatic_model_id = $somatic_build_object->model->genome_model_id;

    my $cn_data_file=$output_dir."/$somatic_model_id.cno_copy_number.csv";
    my $cn_png_file = $cn_data_file.".png";
    
    #run bam-2-cna
    my $job = "gmt somatic bam-to-cna --tumor-bam-file $tumor_bam --normal-bam-file $normal_bam --output-file $cn_data_file --window-size " . $self->window_size;
    my $job_name = $somatic_model_id . "_bam2cna";
    my $oo = $job_name . "_stdout"; #print job's STDOUT in the current directory
    $self->debug_message("Submitting job $job_name (bam-to-cna): $job \n");
    Genome::Sys->shellcmd(cmd => "bsub -q $ENV{GENOME_LSF_QUEUE_BUILD_WORKER} -J $job_name -R 'select[type==LINUX64]' -oo $oo $job");

}

1;
