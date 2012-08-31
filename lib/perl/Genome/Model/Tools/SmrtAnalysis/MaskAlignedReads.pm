package Genome::Model::Tools::SmrtAnalysis::MaskAlignedReads;

use strict;
use warnings;

use Genome;

use File::Temp qw/tempfile/;

class Genome::Model::Tools::SmrtAnalysis::MaskAlignedReads {
    is  => 'Genome::Model::Tools::SmrtAnalysis::Base',
    has_input => [
        cmp_hdf5_file => {
            doc => 'The aligned reads in the format of cmp.h5',
        },
        region_table => {
            doc => 'A FOFN of rgn.h5 tables.',
        },
    ],
    has_optional_input => {
        base_output_directory => {
            is => 'Text',
            doc => 'The base directory to use when running parallel processes.',
        },
        masked_table => {
            doc => 'The output FOFN of masked region tables.',
            is_output => 1,
        },
        log_file => {
            is => 'Text',
            doc => 'Specify a file to log to. Defaults to stderr.',
        },
        debug => {
            is => 'Boolean',
            doc => 'Increases verbosity of logging',
            default_value => 0,
        },
        info => {
            is => 'Boolean',
            doc => 'Display informative log entries',
            default_value => 0,
        },
    },
};


sub help_brief {
    'Mask the aligned reads, typically control sequence.'
}


sub help_detail {
    return <<EOS 
Takes in a rgn.fofn and corresponding cmp.h5. Uses the alignments from the
cmp.h5 to mask corresponding regions of the rgn.h5s. Writes output to a new
rgn.fofn.
EOS
}

sub execute {
    my $self = shift;

    my $cmd = $self->analysis_bin .'/maskAlignedReads.py';
    if ($self->debug) {
        $cmd .= ' --debug';
    }
    if ($self->info) {
        $cmd .= ' --info';
    }
    my @output_files;
    if (defined($self->log_file)) {
        push @output_files, $self->log_file;
        $cmd .= ' --logFile='. $self->log_file;
    }
    if ($self->base_output_directory) {
        if ($self->masked_table) {
            die('Do not define a masked_table when base_output_directory is defined!');
        }
        my (undef, $masked_table) = tempfile('post_control_regions-XXXXXXX', DIR => $self->base_output_directory, SUFFIX => '.fofn');
        $self->masked_table($masked_table);
    }
    push @output_files, $self->masked_table;

    $cmd .= ' '. $self->cmp_hdf5_file .' '. $self->region_table .' '. $self->masked_table;
    $self->shellcmd(
        cmd => $cmd,
        input_files => [$self->cmp_hdf5_file,$self->region_table],
        output_files => \@output_files,
        skip_if_output_is_present => 0,
    );
    return 1;
}


1;
