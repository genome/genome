package Genome::Model::Tools::CopyNumber::IlluminaPipeline;

use strict;
use warnings;
use Genome;
use Cwd;
use FileHandle;

#this is a wrapper to precess illumina file and make copy number graph, swt test. 

class Genome::Model::Tools::CopyNumber::IlluminaPipeline {
    is => 'Command',
    has => [
    output_dir => {
        is => 'String',
        is_optional => 1,
        default => getcwd(),
        doc => 'Directory containing output files and directories.',
    },
    annotation_file => {
        is => 'String',
        is_optional => 1,
        doc => 'Illumina annotation.txt file',
    },
    copy_number_file => {
        is => 'String',
        is_optional => 1,
        doc => 'Illumina copy number output file, typically "...pairedcopynumber.txt".',
    },
    map_file => {
        is => 'String',
        is_optional => 1,
        doc => 'Named map.csv typically. format is [CHR  POS].',
    },
    ]
};

sub help_brief {
    "process Illumina cn files, perform swt test, make cn and swt graphs"
}

sub help_detail {
    "This script processes Illumina copy number data and performs these steps: 1) Creates a map.csv file and per-sample copy number files for a many-sample dataset (unless a map.csv file is already present in the output directory), 2) Creates a copy number plot for each sample, 3) Performs the SWT test for each sample, and 4) creates SWT copy number plots for each sample. SWT plots will be much less noisy than the original cn plots."
}

sub execute {
    my $self = shift;
    my $rlibrary = "cn_lib.R";
    my $output_dir = $self->output_dir;
    my $mapfile = $self->map_file;
    if (-s $mapfile) {
        $self->status_message("A valid map.csv file was input. Skipping the creation of map.csv and per-sample cn files.\n");
    }
    else {
        ## Step 1) Create map.csv file and per-sample copy number file for many-sample datasets
        my $anno_file = $self->annotation_file;
        my $cn_file = $self->copy_number_file;
        unless ($cn_file) {
            $self->error_message("Could not resolve the copy number file from the inputs.");
            return;
        }
        my $create_illumina_cn_files_command = Genome::Model::Tools::CopyNumber::CreateIlluminaCnFiles->create(
            output_dir => $output_dir,
            annotation_file => $anno_file,
            copy_number_file => $cn_file,
        );
        $create_illumina_cn_files_command->execute;
    }

    ## Step 2) Make copy number graph for each sample
    my $cn_graph_command = "per_sample_cn_graph(workdir='$output_dir');";
    my $cn_graph_Rcall = Genome::Model::Tools::R::CallR->create(command=>$cn_graph_command,library=>$rlibrary);
    $cn_graph_Rcall->execute;

    ## Step 3) Perform SWT test
    my $swt_command = "run_swt_test(workdir='$output_dir');";
    my $swt_command_Rcall = Genome::Model::Tools::R::CallR->create(command=>$swt_command,library=>$rlibrary);
    $swt_command_Rcall->execute;

    ##Step 4) Make SWT graphs
    my $swt_graphs_command = "per_sample_swt_graphs(workdir='$output_dir');";
    my $swt_graphs_Rcall = Genome::Model::Tools::R::CallR->create(command=>$swt_graphs_command,library=>$rlibrary);
    $swt_graphs_Rcall->execute;

    return 1;
}
1;
