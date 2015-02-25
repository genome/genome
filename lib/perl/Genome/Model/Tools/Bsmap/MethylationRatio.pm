package Genome::Model::Tools::Bsmap::MethylationRatio;

use strict;
use warnings;
use Cwd;
use Genome;

use Genome::Utility::Text qw(sanitize_string_for_filesystem);
use File::Spec;

my $METHRATIO_COMMAND = File::Spec->join(getcwd, 'methylation_ratio.py');

class Genome::Model::Tools::Bsmap::MethylationRatio {
    is => 'Command',
    has => [
        bam_file => {
            is => 'Text',
            doc => 'An indexed bam file',
            is_input => 1,
        },
        output_file => {
            is => 'Text',
            doc => 'File name for methyl counts (a .bgz extension will be added automatically)',
            is_output => 1,
            is_input => 1,
            default_value => 'methylation',
        },
        threads => {
            is => 'Text',
            doc => 'Numer of threads (2 minimum: one for reading and one for writing; additional threads beyond 2 will be allocated for reading; recommended number is 5)',
            is_output => 0,
            is_input => 1,
            default_value => '2',
        },
        output_directory => {
            is => 'Text',
            doc => 'Where to output the methyl counts',
            is_output => 1,
            is_input => 1,
        },
        reference => {
            is => 'Text',
            doc => '36, 37, or a path to the reference fasta',
            is_input => 1,
        },
    ],
};

sub _reference_fasta {
    
    my ($self) = @_;

    my %mapping = (
        36 => 'NCBI-human-build36',
        37 => 'GRCh37-lite-build37',
    );

    my $reference = $self->reference;

    if (exists $mapping{$reference}) {
        my $reference_build = Genome::Model::Build::ReferenceSequence->get(
            name => $mapping{$reference}
        );
        $reference = $reference_build->cached_full_consensus_path('fa');
    }

    return $reference;
}

sub _generate_command_line {
    my ($self) = @_;

    my @cmd = (
        'python', $METHRATIO_COMMAND,
        '-o', File::Spec->join($self->output_directory, $self->output_file),
        '-d', $self->_reference_fasta,
        '-t', $self->threads
    );
    

    push @cmd, $self->bam_file;

    return join ' ', @cmd;
}

sub execute {
    my ($self) = @_;
    my $return = Genome::Sys->shellcmd(
        cmd => $self->_generate_command_line,
        input_files => [$self->bam_file, $self->_reference_fasta],
    );

    unless ($return) {
        die $self->error_message('Failed to execute: returned %s', $return);
    }
    return 1;
}

1;
