package Genome::Model::Tools::SmrtAnalysis::Pls2Fasta;

use strict;
use warnings;

use Genome;

use File::Temp qw/tempfile/;

class Genome::Model::Tools::SmrtAnalysis::Pls2Fasta {
    is  => 'Genome::Model::Tools::SmrtAnalysis::Base',
    has_input => [
        hdf5_file => {
            doc => 'The pls.h5 or bas.h5 file to convert or an fofn of files',
        },
    ],
    has_optional_input => [
        base_output_directory => {
            is_optional => 1,
            doc => 'This is a base directory where file paths will be resolved for parallel processing only',
        },
        fasta_file => {
            doc => 'The FASTA format file to write output reads',
            is_output => 1,
        },
        region_table => {
            is => 'Text',
            doc => 'The region table used to perform masking or trimming or an fofn of files',
        },
        trim_by_region => {
            is => 'Boolean',
            doc => 'Trim the reads based on the region table.',
            default_value => 1,
        },
        mask_by_region => {
            is => 'Boolean',
            doc => 'Mask the reads based on the region table.',
        },
        split_subreads => {
            is => 'Boolean',
            doc => 'Split the subreads?',
            default_value => 1
        },
        add_sim_coords => {
            is => 'Boolean',
            doc => 'AddSimCoords?',
            default_value => 0,
        },
        min_subread_length => {
            is => 'Number',
            doc => 'The minimum length subread to output in FASTA file.',
        },
    ],
};

sub help_brief {
    'Convert a pls.h5 or bas.h5 file to a FASTA file.'
}


sub help_detail {
    return <<EOS 
    Print reads stored in hdf as fasta.
EOS
}

sub execute {
    my $self = shift;
    my @input_files = ($self->hdf5_file);
    if (defined($self->base_output_directory)) {
        if ($self->fasta_file ) {
            die('Do not specify any output options when base_output_directory is defined!');
        }
        my (undef, $fasta_file) = tempfile('filtered_subreads-XXXXXX', DIR => $self->base_output_directory, SUFFIX => '.fa');
        $self->fasta_file($fasta_file);
    }

    my $cmd = $self->analysis_bin .'/pls2fasta '. $self->hdf5_file .' '. $self->fasta_file;
    if ($self->trim_by_region || $self->mask_by_region) {
        if ($self->trim_by_region && $self->mask_by_region) {
            die('I do not think it would be wise to both trim and mask by region.');
        }
        unless ($self->region_table) {
            die('If you want to trim or mask by region, please define a region table.');
        }
        $cmd .= ' -regionTable '. $self->region_table;
        push @input_files, $self->region_table;
        if ($self->trim_by_region) {
            $cmd .= ' -trimByRegion';
        } elsif ($self->mask_by_region) {
            $cmd .= ' -maskByRegion';
        }
    }
    if (defined($self->min_subread_length)) {
        $cmd .= ' -minSubreadLenght '. $self->min_subread_length;
    }
    unless ($self->split_subreads) {
        $cmd .= ' -noSplitSubreads';
    }
    if ($self->add_sim_coords) {
        $cmd .= ' -addSimCoords';
    }
    $self->shellcmd(
        cmd => $cmd,
        input_files => \@input_files,
        output_files => [$self->fasta_file],
        skip_if_output_is_present => 0,
    );
    return 1;
}


1;
