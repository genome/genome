package Genome::Model::Tools::Bsmap::MethylationRatio;

use strict;
use warnings;
use Cwd;
use Genome;

use Genome::Utility::Text qw(sanitize_string_for_filesystem);
use File::Spec;
use File::Basename qw();

my (undef, $dir) = File::Basename::fileparse(__FILE__);
my $METHRATIO_COMMAND = File::Spec->join($dir, 'methylation_ratio.py');


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
            doc => 'Number of threads (2 minimum: one for reading and one for writing; additional threads beyond 2 will be allocated for reading; recommended number is 5)',
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
		  is => 'Genome::Model::Build::ReferenceSequence',
		  is_input => 1,
		  doc => 'Specify reference build',
		};
    ],
};

sub _reference_fasta {
  my ($self) = @_;
  return $self->reference->cached_full_consensus_path('fa');
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
    Genome::Sys->shellcmd(
        cmd => $self->_generate_command_line,
        input_files => [$self->bam_file, $self->_reference_fasta],
    );

    return 1;
}

1;
