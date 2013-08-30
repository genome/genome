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
        ref_fasta =>  {
            is => 'Path',
        },
    ],
    has_calculated_output => [
        output_dindel_file => {
            is => 'Path',
            calculate => q{ my (undef, undef, $filename) = File::Spec->splitpath($input_vcf);
                            $filename =~ s/.vcf$/.variants.txt/;
                            return File::Spec->join($output_directory, $filename); },
            calculate_from => ['input_vcf', 'output_directory'],
        },
    ],
};

sub help_brief {
    'Convert vcf formatted variants file to dindel format'
}

sub execute {
    my $self = shift;

    $self->create_output_directory();
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
