package Genome::Model::Tools::Sam::Calmd;

use strict;
use warnings;

use Genome;

class Genome::Model::Tools::Sam::Calmd {
    is => ['Genome::Model::Tools::Sam'],
    has => [
        input_file => { },
        output_file => {},
        refseq_file => {},
    ],
};

sub execute {
    my $self = shift;
    my $cmd = $self->samtools_path .' calmd -b '. $self->input_file .' '. $self->refseq_file .' > '. $self->output_file;
    Genome::Sys->shellcmd(
        cmd => $cmd,
        input_files => [$self->input_file,$self->refseq_file],
        output_files => [$self->output_file],
    );
    return 1;
}


1;
