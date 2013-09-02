package Genome::Model::Tools::Dindel::RealignCandidates;

use strict;
use warnings;

use Genome;

class Genome::Model::Tools::Dindel::RealignCandidates {
    is => 'Genome::Model::Tools::Dindel::Base',
    has_input => [
        input_dindel_file => {
            is => 'Path',
        },
        ref_fasta => {
            is => 'Path',
        },
    ],
    has_output => [
        output_dindel_file => {
            is => 'Path',
            is_calculated => 1,
            calculate =>  q{ $output_prefix . ".variants.txt" },
            calculate_from => ['output_prefix'],
        },
    ],
    has_optional_transient => {
        output_prefix => {
            is_calculated => 1,
            calculate =>  q{ File::Spec->join($output_directory, "left_shifted") },
            calculate_from => ['output_directory'],
        },
    },
};

sub help_brief {
    'Left-shift indels in a dindel formatted variants file'
}

sub execute {
    my $self = shift;

    $self->create_output_directory();
    return Genome::Sys->shellcmd_arrayref(
        cmd => [
            $self->dindel_executable,
            '--analysis', 'realignCandidates',
            '--varFile', $self->input_dindel_file,
            '--outputFile', $self->output_prefix,
            '--ref', $self->ref_fasta,
        ],
        input_files => [
            $self->input_dindel_file,
            $self->ref_fasta
        ],
        output_files => [
            $self->output_dindel_file,
        ],
    );
}

1;
