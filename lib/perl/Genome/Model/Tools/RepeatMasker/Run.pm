package Genome::Model::Tools::RepeatMasker::Run;

use strict;
use warnings;

use Genome;

class Genome::Model::Tools::RepeatMasker::Run {
    is => ['Genome::Model::Tools::RepeatMasker','Genome::Sys'],
};

sub execute {
    my $self = shift;

    my $options = ' -species '. $self->species;
    if ($self->mask && $self->mask ne '-n') {
        $options .=  ' '. $self->mask;
    }
    if ($self->sensitivity) {
        $options .= ' '. $self->sensitivity;
    }
    $options .= ' -dir '. $self->output_directory;
    my $cmd = 'RepeatMasker '. $options .' '. $self->fasta_file;
    Genome::Sys->shellcmd(
        cmd => $cmd,
    );
    return 1;
}

1;
