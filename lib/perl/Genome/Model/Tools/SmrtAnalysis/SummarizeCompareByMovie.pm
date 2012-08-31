package Genome::Model::Tools::SmrtAnalysis::SummarizeCompareByMovie;

use strict;
use warnings;

use Genome;

use File::Temp qw/tempfile/;

class Genome::Model::Tools::SmrtAnalysis::SummarizeCompareByMovie {
    is  => 'Genome::Model::Tools::SmrtAnalysis::Base',
    has_input => [
        cmp_hdf5_file => {
            doc => 'The aligned reads in the format of cmp.h5',
        },
    ],
    has_optional_input => {
        base_output_directory => {
            is => 'Text',
            doc => 'The base output directory used for parallel processing',
        },
        output_csv_file => {
            doc => 'The output summary report in csv format.',
            is_output => 1,
        },
        original => {
            doc => 'Original reads FASTA file',
        },
        filtered => {
            doc => 'Filtered reads FASTA file',
        },
        fastq => {
            doc => 'Filtered reads FASTQ file',
        },
        fofn => {
            doc => 'File of bas.h5 files',
        },
        external => {
            is => 'Boolean',
            doc => 'Only output PacBio-external metrics',
            default_value => 1,
        },
    },
};


sub help_brief {
    'Breaks out post-mapping statistics by movie.'
}


sub help_detail {
    return <<EOS
Breaks out post-mapping statistics by movie.
EOS
}

sub execute {
    my $self = shift;

    my $cmd = $self->analysis_bin .'/summarizeCompareByMovie.py';
    my @input_files;
    if (defined($self->original)) {
        $cmd .= ' --original='. $self->original;
        push @input_files, $self->original;
    }
    if (defined($self->filtered)) {
        $cmd .= ' --filtered='. $self->filtered;
        push @input_files, $self->filtered;
    }
    if (defined($self->fastq)) {
        $cmd .= ' --fastq='. $self->fastq;
        push @input_files, $self->fastq;
    }
    if (defined($self->fofn)) {
        $cmd .= ' --fofn='. $self->fofn;
        push @input_files, $self->fofn;
    }
    if ($self->external) {
        $cmd .= ' --external';
    }
    if ($self->base_output_directory) {
        if ($self->output_csv_file) {
            die('Do not define an output_csv_file when base_output_directory is defined!');
        }
        my (undef, $output_csv_file) = tempfile('control_results_by_movie-XXXXXXX', DIR => $self->base_output_directory, SUFFIX => '.csv');
        $self->output_csv_file($output_csv_file);
    }
    $cmd .= ' '. $self->cmp_hdf5_file .' > '. $self->output_csv_file;
    $self->shellcmd(
        cmd => $cmd,
        input_files => \@input_files,
        output_files => [$self->output_csv_file],
        skip_if_output_is_present => 0,
    );
    return 1;
}


1;
