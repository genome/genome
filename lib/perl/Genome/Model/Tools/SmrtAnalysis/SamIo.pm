package Genome::Model::Tools::SmrtAnalysis::SamIo;

use strict;
use warnings;

use Genome;
use File::Basename;

class Genome::Model::Tools::SmrtAnalysis::SamIo {
    is => ['Genome::Model::Tools::SmrtAnalysis::Base'],
    has_input => {
        cmp_hdf5_file => {
            doc => 'The input cmp.h5 file.',
            is_output => 1,
        },
    },
    has_optional_input => {
        output_file => {
            is => 'Text',
            doc => 'Name of output SAM file',
            is_output => 1,
        },
        reference_repository => {
            is => 'Text',
            doc => 'Path to Reference Sequence Repository',
        },
        info => {
            is => 'Boolean',
            doc => 'Turn on progress monitoring (stdout)',
            default_value => 0,
        },
        bam => {
            is => 'Boolean',
            doc => 'Generate BAM file',
            default_value => 1,
        },
        md => {
            is => 'Boolean',
            doc => 'Populate MD field in generated SAM file',
            default_value => 1,
        },
    },
    has_optional_param => [
        lsf_resource => { default_value => "-g /pacbio/smrtanalysis -M 16000000 -R 'select[type==LINUX64 && mem>=16000 && tmp>=160000] rusage[mem=16000,tmp=80000]'" },
    ],
};

sub help_brief {
    'Convert a cmp.h5 file to a SAM/BAM file.'
}

sub help_detail {
    return <<EOS
Convert a cmp.h5 file to a SAM/BAM file.
EOS
}

sub execute {
    my $self = shift;

    my $cmd = $self->analysis_bin .'/SAMIO.py';
    if ($self->info) {
        $cmd .= ' --info';
    }
    if ($self->bam) {
        $cmd .= ' --bam';
    }
    if ($self->md) {
        $cmd .= ' --md';
    }
    my @input_files;
    my @output_files;
    if (defined($self->reference_repository)) {
        $cmd .= ' --refrepos='. $self->reference_repository;
        push @input_files, $self->reference_repository;
    }
    unless (defined($self->output_file)) {
        my ($basename, $dirname, $suffix) = File::Basename::fileparse($self->cmp_hdf5_file,qw/\.cpm\.h5/);
        my $output_file = $dirname .'/'. $basename .'.sam';
        $self->output_file($output_file);
    }
    $cmd .= ' --outfile='. $self->output_file;
    push @output_files, $self->output_file;

    $cmd .= ' '. $self->cmp_hdf5_file;
    push @input_files, $self->cmp_hdf5_file;
    $self->shellcmd(
        cmd => $cmd,
        input_files => \@input_files,
        output_files => \@output_files,
        skip_if_output_is_present => 0,
    );
    return 1;
};


1;
