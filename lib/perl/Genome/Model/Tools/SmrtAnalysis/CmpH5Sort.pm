package Genome::Model::Tools::SmrtAnalysis::CmpH5Sort;

use strict;
use warnings;

use Genome;

class Genome::Model::Tools::SmrtAnalysis::CmpH5Sort {
    is => ['Genome::Model::Tools::SmrtAnalysis::Base'],
    has_input => {
        input_cmp_hdf5_file => {
            doc => 'The input cmp.h5 file that can also be the output if a output file is not defined',
            is_output => 1,
        },
    },
    has_optional_input => {
        output_cmp_hdf5_file => {
            doc => 'The output cmp hdf5 file that has been sorted',
            is_output => 1,
        },
        silent => {
            is => 'Boolean',
            doc => 'print nothing.',
            default_value => 0,
        },
        verbose => {
            is => 'Boolean',
            doc => 'print debugging information.',
            default_value => 0,
        },
        deep => {
            is => 'Boolean',
            doc => 'whether a deep sorting should be conducted, i.e. sort the AlignmentArrays',
            default_value => 0,
        },
        jobs => {
            doc => 'Number of child processes to launch. This only speeds up processing if there are multiple references groups. Not yet Implemented.',
        },
        tmp_dir => {
            doc => 'Temporary directory to use when sorting in-place.',
        },
    },
};

sub help_brief {
    'Sort cmp.h5 files.'
}

sub help_detail {
    return <<EOS
Sort cmp.h5 files. If output-file is unspecified the input-file is
overwritten. If there are a number of reference groups then the
indexing processing can occur in parallel.
EOS
}

sub execute {
    my $self = shift;

    my $cmd = $self->analysis_bin .'/cmpH5Sort.py';
    if ($self->silent) {
        $cmd .= ' --silent';
    }
    if ($self->verbose) {
        $cmd .= ' --verbose';
    }
    if ($self->deep) {
        $cmd .= ' --deep';
    }
    if (defined($self->jobs)) {
        $cmd .= ' --jobs='. $self->jobs;
    }
    if (defined($self->tmp_dir)) {
        $cmd .= ' --tmpDir='. $self->tmp_dir;
    } else {
        my $tmp_dir = File::Temp::tempdir(
            'Genome-Model-Tools-CmpH5Sort-XXXXXX',
            DIR => '/gscmnt/gc2123/production/lsf_shared_dir',
            CLEANUP => 1,
        );;
        $cmd .= ' --tmpDir='. $tmp_dir;
    }
    $cmd .= ' '. $self->input_cmp_hdf5_file;
    my @output_files;
    if (defined($self->output_cmp_hdf5_file)) {
        $cmd .= ' '. $self->output_cmp_hdf5_file;
        push @output_files, $self->output_cmp_hdf5_file;
    }
    $self->shellcmd(
        cmd => $cmd,
        input_files => [$self->input_cmp_hdf5_file],
        output_files => \@output_files,
        skip_if_output_is_present => 0,
    );
    return 1;
}


1;
