package Genome::Model::Tools::Dindel::GetCigarIndels;

use strict;
use warnings;

use Genome;

class Genome::Model::Tools::Dindel::GetCigarIndels {
    is => 'Genome::Model::Tools::Dindel::Base',
    has_input => [
        input_bam => {
            is => 'Path',
            doc => 'The .bai file must be right next to it.',
        },
        ref_fasta => {
            is => 'Path',
        },
    ],
    has_output => [
        output_variants => {
            is => 'Path',
            is_calculated => 1,
            calculate =>  q{ "$output_prefix.variants.txt" },
            calculate_from => ['output_prefix'],
        },
        output_libraries => {
            is => 'Path',
            is_calculated => 1,
            calculate =>  q{ "$output_prefix.libraries.txt" },
            calculate_from => ['output_prefix'],
        },
    ],
    has_optional_transient => {
        output_prefix => {
            is_calculated => 1,
            calculate =>  q{ File::Spec->join($output_directory, "cigar_generated_indels") },
            calculate_from => ['output_directory'],
        },
    },
};

sub help_brief {
    'Call indels in the input_bam (output_variants) and calculate the insert-size-distribution (output_libraries)'
}

sub execute {
    my $self = shift;

    $self->create_output_directory();
    return Genome::Sys->shellcmd_arrayref(
        cmd => [
            $self->dindel_executable,
            '--analysis', 'getCIGARindels',
            '--bamFile', $self->input_bam,
            '--ref', $self->ref_fasta,
            '--outputFile', $self->output_prefix,
        ],
        input_files => [
            $self->input_bam,
            $self->ref_fasta,
        ],
        output_files => [
            $self->output_variants,
            $self->output_libraries,
        ],
    );
}

1;
