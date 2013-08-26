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
        output_prefix => {
            is => 'String',
        },
    ],
};

sub help_brief {
    'Run getCIGARindels'
}

sub help_synopsis {
    return <<EOS
EOS
}

sub help_detail {
    return <<EOS
EOS
}


sub execute {
    my $self = shift;

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
    );
}

1;
