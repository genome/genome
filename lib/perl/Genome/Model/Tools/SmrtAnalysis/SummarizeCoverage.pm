package Genome::Model::Tools::SmrtAnalysis::SummarizeCoverage;

use strict;
use warnings;

use Genome;
use File::Basename;

class Genome::Model::Tools::SmrtAnalysis::SummarizeCoverage {
    is => ['Genome::Model::Tools::SmrtAnalysis::Base'],
    has_input => {
        cmp_hdf5_file => {
            doc => 'The input cmp.h5 file.',
            is_output => 1,
        },
        summary_file => {
            is => 'Text',
            doc => 'The output GFF3 format summary file.',
            is_output => 1,
        },
    },
    has_optional_input => {
        reference_repository => {
            is => 'Text',
            doc => 'Reference Repository for mapping fullName --> RepoID',
        },
        number_of_regions => {
            is => 'Number',
            doc => 'Desired number of regions in the summary statistics(used for guidance, not strict)',
            default_value => 500,
        },
        region_size => {
            is => 'Number',
            doc => 'If supplied used a fixed region size [overrides numRegions]',
        },
        debug => {
            is => 'Boolean',
            doc => 'Turn on debugging output',
            default_value => 0,
        },
    },
};

sub help_brief {
    'Generate a GFF3 file summarizing various aspects of multi-read alignment(cmp.h5).'
}

sub help_detail {
    return <<EOS
Generate a GFF3 file summarizing various aspects of multi-read alignment(cmp.h5).
EOS
}

sub execute {
    my $self = shift;

    my $cmd = $self->analysis_bin .'/summarizeCoverage.py';
    if ($self->debug) {
        $cmd .= ' --debug';
    }
    if (defined($self->region_size)) {
        $cmd .= ' --regionSize='. $self->region_size;
    }
    if (defined($self->number_of_regions)) {
        $cmd .= ' --numRegions='. $self->number_of_regions;
    }
    my @input_files;
    my @output_files;
    if (defined($self->reference_repository)) {
        $cmd .= ' --reference='. $self->reference_repository;
        push @input_files, $self->reference_repository;
    }

    push @output_files, $self->summary_file;

    $cmd .= ' '. $self->cmp_hdf5_file .' > '. $self->summary_file;
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
