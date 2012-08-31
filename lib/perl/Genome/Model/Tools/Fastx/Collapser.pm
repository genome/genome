package Genome::Model::Tools::Fastx::Collapser;

use strict;
use warnings;

use Genome;
use Genome::Sys;
use File::Basename;

class Genome::Model::Tools::Fastx::Collapser {
    is => ['Genome::Model::Tools::Fastx'],
    has_constant => {
        fastx_tool => { value => 'fastx_collapser' },
    },
    has_input => [
        input_file => {
            is => 'Text',
            doc => 'The input FASTQ/A file to collapse.(This works on fastq but I get errors about the quality string)',
        },
    ],
    has_output => [
        output_file => {
            is => 'Text',
            doc => 'The output FASTQ/A file containing collapsed sequence.',
            is_optional => 1,
        }
    ],
};

sub execute {
    my $self = shift;
    unless (Genome::Sys->validate_file_for_reading($self->input_file)) {
        $self->error_message('Failed to validate fastq file for read access '. $self->input_file .":  $!");
        die($self->error_message);
    }
    my @suffix = qw/fq fa fastq fasta fna txt/;
    my ($basename,$dirname,$suffix) = File::Basename::fileparse($self->input_file,@suffix);
    $basename =~ s/\.$//;
    unless ($self->output_file) {
        $self->output_file($dirname .'/'. $basename .'_collapsed.fa');
    }
    unless (Genome::Sys->validate_file_for_writing($self->output_file)) {
        $self->error_message('Failed to validate output file for write access '. $self->output_file .":  $!");
        die($self->error_message);
    }
    #fastx_quality_stats (as of 0.0.7) won't process from tmp with -i and -o
    my $cmd = $self->fastx_tool_path .' -i '. $self->input_file .' -o '. $self->output_file;
    Genome::Sys->shellcmd(
        cmd => $cmd,
        input_files => [$self->input_file],
        output_files => [$self->output_file],
    );
    return 1;
}

1;
