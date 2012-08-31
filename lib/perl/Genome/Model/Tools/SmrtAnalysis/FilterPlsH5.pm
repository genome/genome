package Genome::Model::Tools::SmrtAnalysis::FilterPlsH5;

use strict;
use warnings;

use Genome;
use File::Temp qw/tempfile/;


class Genome::Model::Tools::SmrtAnalysis::FilterPlsH5 {
    is  => 'Genome::Model::Tools::SmrtAnalysis::Base',
    has_input => [
        input_fofn => {
            doc => 'Specify the location of a FOFN file pointing to the input h5 files to be filtered.',
            is_output => 1,
        },
    ],
    has_optional_input => [
        base_output_directory => {
            is_optional => 1,
            doc => 'This is a base directory where file paths will be resolved for parallel processing only',
        },
        output_summary => {
            is_output => 1,
            doc => 'Optional name of a file to write filtering information to.',
        },
        output_fofn => {
            is_output => 1,
            doc => 'Specify the location of a FOFN file pointing to the regions files.',
        },
        output_directory => {
            is_output => 1,
            doc => 'Specify a directory to write results to in the format of regions tables (rgn.h5s).',
        },
        debug => {
            doc => 'Outputs a log to stderr with helpful debug info.',
        },
        log_file => {
            is_output => 1,
            doc => 'Set a log file for logging output. Defaults to stderr.',
        },
        trim => {
            is => 'Boolean',
            default_value => 1,
            doc => 'Trim down to HQRegion.',
        },
        max_sandwiches => {
            is => 'Number',
            doc => 'Maximum # of Sandwiches allowed(0.0-Inf)',
        },
        min_read_score => {
            is => 'Number',
            doc => 'Accuracy Prediction from Primary(0.0-1.0)',
        },
        min_read_length => {
            is => 'Number',
            doc => 'Minimum # of Pulses in the ZMW(0.0-Inf)',
        },
        read_white_list => {
            is => 'Text',
            doc => 'File path to newline separated list of ReadIds',
        },
        min_snr => {
            is => 'Number',
            doc => 'Minimum SNR cutoff per Read(0.0-Inf)',
        },
    ],
};

sub help_brief {
    "Filter a pls.h5 or bas.h5 file down to high-quality regions.",
}


sub help_detail {
    return <<EOS 
This script takes in a pls.h5 FOFN (.fofn) plus one or more filter
specifications. A filter specification takes the form: filterName:threshold.
For a list of filters and valid thresholds, use the --availableFilters option.
Also note that trimming by HQRegions is on by default, and can be disabled.
EOS
}


sub execute {
    my $self = shift;

    my @filters;
    if (defined($self->max_sandwiches)) {
        push @filters, 'MaxSandwiches='. $self->max_sandwiches;
    }
    if (defined($self->min_read_score)) {
        push @filters, 'MinReadScore='. $self->min_read_score;
    }
    if (defined($self->min_read_length)) {
        push @filters, 'MinRL='. $self->min_read_length;
    }
    if (defined($self->read_white_list)) {
        push @filters, 'ReadWhitelist='. $self->read_white_list;
    }
    if (defined($self->min_snr)) {
        push @filters, 'MinSNR='. $self->min_snr;
    }
    unless (@filters) {
        die('Failed to define one filter!');
    }
    my $filters = join(',',@filters);
    my $cmd = $self->analysis_bin ."/filterPlsH5.py --filter='$filters'";
    if ($self->trim) {
        $cmd .= ' --trim=True';
    } else {
        $cmd .= ' --trim=False';
    }
    if (defined($self->base_output_directory)) {
        if ($self->output_directory || $self->output_fofn || $self->output_summary) {
            die('Do not specify any output options when base_output_directory is defined!');
        }
        $self->output_directory($self->base_output_directory .'/filtered_regions');
        my (undef, $output_summary) = tempfile('filtered_summary-XXXXXX', DIR => $self->base_output_directory, SUFFIX => '.csv');
        $self->output_summary($output_summary);
        my (undef, $output_fofn) = tempfile('filtered_regions-XXXXXX', DIR => $self->base_output_directory, SUFFIX => '.fofn');
        $self->output_fofn($output_fofn);
    }
    if (defined($self->output_directory)) {
        unless (-d $self->output_directory) {
            unless (Genome::Sys->create_directory($self->output_directory)) {
                die('Failed to created directory: '. $self->output_directory);
            }
        }
        $cmd .= ' --outputDir='. $self->output_directory;
    }
    my @output_files;
    if (defined($self->output_fofn)) {
        push @output_files, $self->output_fofn;
        $cmd .= ' --outputFofn='. $self->output_fofn;
    }
    if (defined($self->output_summary)) {
        push @output_files, $self->output_summary;
        $cmd .= ' --outputSummary='. $self->output_summary;
    }
    if ($self->debug) {
        $cmd .= ' --debug';
    }
    if (defined($self->log_file)) {
        push @output_files, $self->log_file;
        $cmd .= ' --logFile='. $self->log_file;
    }

    my $input_fofn_fh = Genome::Sys->open_file_for_reading($self->input_fofn);
    for my $line ($input_fofn_fh) {
        chomp($line);
        if ($line =~ /^\s*$/) { next; }
        my $file = $line;
        unless (-f $file && -s $file) {
            die('Failed to find input file: '. $file );
        }
    }
    $input_fofn_fh->close;
    $cmd .= ' '. $self->input_fofn;

    # TODO: Resolve the input and output files/directories
    $self->shellcmd(
        cmd => $cmd,
        input_files => [$self->input_fofn],
        output_files => \@output_files,
        skip_if_output_is_present => 0,
    );

    return 1;
}


1;
