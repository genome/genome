package Genome::Model::Tools::Dindel::RealignCandidates;

use strict;
use warnings;

use Genome;

class Genome::Model::Tools::Dindel::RealignCandidates {
    is => 'Command',
    has_input => [
        variant_file => {
            is => 'Path',
        },
        output_prefix => {
            is => 'String',
        },
        ref_fasta => {
            is => 'Path',
        },
    ],
    has_output => [
        output_file => {
            is => 'Path',
            is_calculated => 1,
            calculate =>  q{ "$output_prefix.variants.txt" },
            calculate_from => ['output_prefix'],
        }
    ]
};

sub help_brief {
    'Left-shift indels you got from a vcf'
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
    my $dindel_location = "/gscmnt/gc2146/info/medseq/dindel/binaries/dindel-1.01-linux-64bit";
    return Genome::Sys->shellcmd_arrayref(
        cmd => [
            $dindel_location,
            '--analysis', 'realignCandidates',
            '--varFile', $self->variant_file,
            '--outputFile', $self->output_prefix,
            '--ref', $self->ref_fasta,
        ],
        input_files => [
            $self->variant_file,
            $self->ref_fasta
        ],
    );
}

1;
