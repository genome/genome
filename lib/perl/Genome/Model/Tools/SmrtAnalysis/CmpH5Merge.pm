package Genome::Model::Tools::SmrtAnalysis::CmpH5Merge;

use strict;
use warnings;

use Genome;

class Genome::Model::Tools::SmrtAnalysis::CmpH5Merge {
    is => ['Genome::Model::Tools::SmrtAnalysis::Base'],
    has_input => {
        input_cmp_hdf5_files => {
            doc => 'The input cmp.h5 file that can also be the output if a output file is not defined',
        },
        output_cmp_hdf5_file => {
            doc => 'The output cmp hdf5 file that has been sorted',
            is_output => 1,
        },
    },
    has_optional_input => {
        force_merge => {
            is => 'Boolean',
            doc => 'Bypass validation of cmp.h5 files before merging and force a merge',
            default_value => 0,
        },
        info => {
            is => 'Boolean',
            doc => 'Display informative log entries',
            default_value => 0,
        },
        debug => {
            is => 'Boolean',
            doc => 'Increases verbosity of logging',
            default_value => 0,
        },
        log_file => {
            doc => 'Specify a file to log to. Defaults to stderr.',
        },
        remove_originals => {
            default_value => 1,
        },
    },
};

sub help_brief {
    'Merge cmp.h5 files.'
}

sub help_detail {
    return <<EOS
Merges one or more CmpH5 files into a single, equivalent CmpH5 file.
EOS
}

sub execute {
    my $self = shift;

    my $input_files = $self->input_cmp_hdf5_files;
    unless (ref($input_files) eq 'ARRAY') {
        my @input_files = split(',',$input_files);
        unless (scalar(@input_files)) {
            die('Must define one or more input_cmp_hdf5_files as a comma separated list!');
        }
        $input_files = \@input_files;
    }
    unless (scalar(@{$input_files})) {
        die('Must define one or more input_cmp_hdf5_files!');
    }
    my $cmd = $self->analysis_bin .'/cmpH5Merge.py';
    if ($self->debug) {
        $cmd .= ' --debug';
    }
    if ($self->info) {
        $cmd .= ' --info';
    }
    if ($self->force_merge) {
        $cmd .= ' --forceMerge';
    }
    my @output_files;
    if (defined($self->log_file)) {
        $cmd .= ' --logFile='. $self->log_file;
        push @output_files, $self->log_file;
    }
    $cmd .= ' '. join(' ', @{$input_files}) .' '. $self->output_cmp_hdf5_file;
    push @output_files, $self->output_cmp_hdf5_file;
    print $cmd ."\n";
    $self->shellcmd(
        cmd => $cmd,
        input_files => $input_files,
        output_files => \@output_files,
        skip_if_output_is_present => 0,
    );
    if ($self->remove_originals) {
        for my $input_file (@{$input_files}) {
            unlink($input_file) || die ('Failed to remove original file '. $input_file);
        }
    }
    return 1;
}


1;
