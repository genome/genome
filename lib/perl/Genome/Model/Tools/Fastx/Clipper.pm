package Genome::Model::Tools::Fastx::Clipper;

use strict;
use warnings;

use Genome;
use Genome::Sys;
use File::Basename;

class Genome::Model::Tools::Fastx::Clipper {
    is => ['Genome::Model::Tools::Fastx'],
    has_constant => {
        fastx_tool => { value => 'fastx_clipper' },
    },
    has_input => [
        input_file => {
            doc => 'The input FASTQ/A file to collapse.(This works on fastq but I get errors about the quality string)',
        },
        output_directory => {
            is_optional => 1,
        },
        params => { },
    ],
    has_output => [
        output_file => {
            doc => 'The output FASTQ/A file containing collapsed sequence.',
            is_optional => 1,
        },
        log_file => {
            is_optional => 1,
        },
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
        $self->output_file($dirname .'/'. $basename .'_clipped.'. $suffix);
    }
    unless ($self->log_file) {
        $self->log_file($dirname .'/'. $basename .'_clipped.log');
    }
    unless (Genome::Sys->validate_file_for_writing($self->output_file)) {
        $self->error_message('Failed to validate output file for write access '. $self->output_file .":  $!");
        die($self->error_message);
    }
    unless (Genome::Sys->validate_file_for_writing($self->log_file)) {
        $self->error_message('Failed to validate output file for write access '. $self->log_file .":  $!");
        die($self->error_message);
    } 
    my $params = $self->params . ' -v ';
    my $cmd = $self->fastx_tool_path .' '. $params .' -i '. $self->input_file .' -o '. $self->output_file .' > '. $self->log_file;
    Genome::Sys->shellcmd(
        cmd => $cmd,
        input_files => [$self->input_file],
        output_files => [$self->output_file,$self->log_file],
    );
    return 1;
}

1;
