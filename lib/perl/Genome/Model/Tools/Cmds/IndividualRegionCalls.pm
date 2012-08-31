package Genome::Model::Tools::Cmds::IndividualRegionCalls;

use warnings;
use strict;
use Genome;
use Cwd;
use Statistics::R;
require Genome::Sys;

class Genome::Model::Tools::Cmds::IndividualRegionCalls {
    is => 'Command',
    has => [
    cmds_test_dir => {
        type => 'String',
        is_optional => 0,
        is_input => 1,
        doc => 'Directory containing CMDS .test result files (named cmds_test by gmt cmds execute).',
    },
    cmds_input_data_dir => {
        type => 'String',
        is_optional => 0,
        is_input => 1,
        doc => 'Directory containing original input data files (only!) to gmt cmds execute. Files must be labeled as filename.chr (convention of other CMDS tools).',
    },
    output_dir => {
        type => 'String',
        is_optional => 1,
        is_input => 1,
        is_output => 1,
        default => getcwd(),
        doc => 'Directory for output of ROI file and individual call files showing per-sample variations. Default: current working directory.',        
    },
    cmds_test_cutoff => {
        type => 'Number',
        is_optional => 1,
        default => 3,
        doc => 'CMDS test cutoff value (-log10(z.p in cmds.test output))',
    },
    window_size => {
        type => 'Number',
        is_optional => 1,
        default => 100000,
        doc => 'Look for p-values passing the cutoff and group all values that are within this far of each other.',
    },
    permutations => {
        type => 'Number',
        is_optional => 1,
        default => "NA",
        doc => 'Optional number of permutations to give weight to the per-sample test results of ROI CNA alteration. More permutations give more confidence in calls, but will take more time to produce.',
    },
    ]
};

sub help_brief {
    'Look at all individuals in the study set at significant cmds_test regions'
}

sub help_detail {
    "This script takes the CMDS results and looks for regions of interest (ROIs) which show significant copy-number alterations. Significance is based on the CMDS test value, -log10(p-value), where the default cutoff = 3. Windows with significant p-values are searched and grouped if within 100,000 bp (default). These ROIs are saved in the output_dir in a file called ROIs.txt. This ROI file is then parsed to create individual region calls using Qunyuan's regioncall.R script. regioncall.R takes each ROI and creates an output file with per-sample CNA information in this region. THese output files are put in the --output-dir."
}

sub execute {
    my $self = shift;
    my $cmds_test_dir = $self->cmds_test_dir;
    my $data_dir = $self->cmds_input_data_dir;
    my $output_dir = $self->output_dir;
    my $permutations = $self->permutations;
    my $cmds_test_cutoff = $self->cmds_test_cutoff;
    my $pv_cutoff = exp(-log(10)*$cmds_test_cutoff); #used to compare with the z.p column of cmds.test output
    my $window_frame = $self->window_size;
    my %regions; # to store data to be used with R calls below
    my @chrs; # to store chrs for R calls

    #open output file
    my $outfile = $output_dir . "/ROIs.txt";
    my $out_fh = new IO::File $outfile,"w";

    #loop through all cmds output .test files
    opendir CMDS_TEST,$cmds_test_dir;
    while (my $file = readdir CMDS_TEST) {
        next if ($file eq "." || $file eq "..");
        my $full_path_file = "$cmds_test_dir/$file";
        my $cmds_test_fh = new IO::File $full_path_file,"r" or die "Cannot open file $file\n";
        (my $chr = $file) =~ s/.+\.(.+)\.test/$1/;
        my $cur_reg_start; # used to record expanding window as more data meeting cutoff is found
        my $cur_reg_end; # used to record expanding window as more data meeting cutoff is found

        #parse each line of each file
        while (my $line = $cmds_test_fh->getline) {
            next if ($line =~ /^chromosome/);
            my ($chr,$window,$start,$mid,$end,$m,$z,$m_sd,$m_p,$z_sd,$z_p) = split /\s+/, $line;
            next if ($z_p eq "NA");
            if ($z_p <= $pv_cutoff) {
                if (!defined $cur_reg_start) {
                    $cur_reg_start = $start;
                    $cur_reg_end = $end;
                    next;
                }
                my $distance = $end - $cur_reg_end;
                if ($distance <= $window_frame) {
                    $cur_reg_end = $end;
                    next;
                }
                if ($distance > $window_frame) {
                    $regions{$chr}{$cur_reg_start} = $cur_reg_end;
                    $cur_reg_start = $start;
                    $cur_reg_end = $end;
                }
            }#end, z.p passes cutoff
        }#end, parsing cmds_test file

        #capture last region left in queue if there was one
        if (defined $cur_reg_start) {
            $regions{$chr}{$cur_reg_start} = $cur_reg_end;
        }

        #write results to ROIs.txt file (if there were significant regions)
        @chrs = keys %regions;
        if (scalar @chrs) {
            my @region_starts = keys %{ $regions{$chr} };
            for my $start (@region_starts) {
                $out_fh->print("$chr\t$start\t$regions{$chr}{$start}\t$permutations\n");
            }
        }
    }#end, chr files
    $out_fh->close;

    #read cmds input data dir 
    opendir DATA,$data_dir;
    while (my $file = readdir DATA ) {
        next if ($file eq "." || $file eq "..");
        (my $chr = $file) =~ s/.+\.(.+)$/$1/;

        #if $chr from data file matches a signficant region, perform Region_Call
        if (grep { /^$chr$/ } @chrs) {
            my $full_path_data_file = $data_dir . "/" . $file;
            for my $start (keys %{ $regions{$chr} }) {
                my $command = "Region_calls(datafile='$full_path_data_file',chr='$chr',start=$start,end='$regions{$chr}{$start}',permun=$permutations,output_dir='$output_dir');";
                my $result_object = Genome::Model::Tools::R::CallR->execute(command => $command, library => "cmds_lib.R");
                unless($result_object->result == 1){
                    $self->error_message("callR with $command failed - result = " . $result_object->result);
                    die $self->error_message;
                }
            }
        }
    }#end, reading original data dir

    return 1;
}
1;
