package Genome::Model::Tools::Cmds::Execute;

use warnings;
use strict;
use Genome;
use Cwd;

class Genome::Model::Tools::Cmds::Execute {
    is => 'Command',
    has => [
    data_directory => {
        type => 'String',
        is_optional => 0,
        is_input => 1,
        is_output => 1,
        doc => 'Directory containing all of the input files for sample group (and nothing else!!).',
    },
    window_size => {
        type => 'String',
        is_optional => 0,
        is_input => 1,
        default => 30,
        doc => 'Window size (# of markers) for CMDS computations.',
    },
    step_size => {
        type => 'String',
        is_optional => 0,
        is_input => 1,
        default => 1,
        doc => 'Window step size for CMDS computations.',
    },
    output_directory => {
        type => 'String',
        is_optional => 1,
        is_input => 1,
        default => getcwd(),
        doc => 'Directory for output folders cmds_test and cmds_plot to be created. Default: current working directory.',
        
    },
    test_output_directory => {
        type => 'String',
        is_output => 1,
        calculate_from => ['output_directory'],
        calculate => q{ return $output_directory . "/cmds_test"; },
        doc => 'Directory containing the .test output files',
    },
    plot_output_directory => {
        type => 'String',
        calculate_from => ['output_directory'],
        calculate => q{ return $output_directory . "/cmds_plot"; },
        doc => 'Directory containing the plot output files',
    },
    ],
};

sub help_brief {
    "Run CMDS analysis"
}

sub help_detail {
    "This script runs the CMDS tool written by Qunyuan. It submits a job array to parallelize the process by chromosome. The process runs on every file in the data directory, so make sure it is clean. It creates two folders for output in the 'output_directory', one called cmds_test, containing text files of results, and a second caled cmds_plot, containing results in the form of images. DO NOT RUN THIS SCRIPT IN THE DATA DIRECTORY."
}

sub execute {
    my $self = shift;

    unless (`uname -a` =~ /x86_64/){
        $self->error_message("This tool produces different results when run on a 32 bit system. It can be done, but results will differ from running on a 64 bit system. Please run on 64-bit for the sake of consistency.");
        die;
    }

    #define input and output directories
    my $data_dir = $self->data_directory;
    my $output_dir = $self->output_directory;

    #gmt r call-r will create both cmds_test and cmds_plot directory for us
    my $plot_dir = $self->plot_output_directory;
    my $test_dir = $self->test_output_directory;
    
    #submit an R job for each chromosome input file
    my $index = 1;
    opendir DATA_DIR,$data_dir;
    while (my $file = readdir DATA_DIR) {
        next if ($file eq "." || $file eq "..");
        my $command = "cmds.focal.test(data.dir='$data_dir',wsize=" . $self->window_size . ",wstep=" . $self->step_size . ",analysis.ID='$index',chr.colname='CHR',pos.colname='POS',plot.dir='$plot_dir',result.dir='$test_dir');";
        my $callr = Genome::Model::Tools::R::CallR->create(
            command => $command,
            library => 'cmds_lib.R',
        );
        my $result = $callr->execute;
        unless ($result == 1) {
            $self->error_message("CallR failed for command: $command");
            die;
        }
        $index++;
    }

    return 1;
}
1;
