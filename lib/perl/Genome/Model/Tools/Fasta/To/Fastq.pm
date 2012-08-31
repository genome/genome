package Genome::Model::Tools::Fasta::To::Fastq;

use strict;
use warnings;

use Genome;

class Genome::Model::Tools::Fasta::To::Fastq {
    is           => 'Genome::Model::Tools::Fasta::To',
    has_optional => [
        dir  => {
            is      => 'String',
            doc     => 'The path to place fastq file',
            default => '.',
        },
        out_filename => {
            is      => 'String',
            doc     => 'The name of output fastq file',
        },
    ],
};


sub help_brief {
    "convert fasta sequence (and qual files) to fastq files"
}


sub help_detail {                           
    return <<EOS 
This tool convert fasta sequence (and quality values) to fastq files using BioPerl.
EOS
}


sub out_file {
    my $self = shift;
    return $self->dir.'/'.$self->out_filename if $self->out_filename;
    return $self->dir.'/'.$self->fasta_basename.'.'.$self->_format_type;
}


sub _format_type {
    return 'fastq';
}


sub write_method {
    return 'write_fastq';
}




1;
