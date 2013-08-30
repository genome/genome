package Genome::Model::Tools::Dindel::VcfToDindel;

use strict;
use warnings;

use Genome;

class Genome::Model::Tools::Dindel::VcfToDindel {
    is => 'Genome::Model::Tools::Dindel::Base',
    has_input => [
        input_vcf =>  {
            is => 'Path',
        },
        output_dindel_file => {
            is => 'Path',
            is_output => 1,
        },
        ref_fasta =>  {
            is => 'Path',
        },
    ],
};

sub help_brief {
    'Turn vcf into stupid dindel format'
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

    my @cmd = (
        'python', $self->python_script('convertVCFToDindel'),
        '--inputFile', $self->input_vcf,
        '--outputFile', $self->output_dindel_file,
        '--refFile', $self->ref_fasta,
    );

    return Genome::Sys->shellcmd_arrayref(
        cmd => \@cmd,
        input_files => [
            $self->input_vcf,
            $self->ref_fasta,
        ],
        output_files => [
            $self->output_dindel_file,
        ],
    );
}

1;
