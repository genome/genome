package Genome::Model::Tools::Fasta::Dust;

use strict;
use warnings;

use Genome;
use File::Basename;
use IO::File;

class Genome::Model::Tools::Fasta::Dust {
    is => 'Genome::Model::Tools::Fasta',
    has_input => [
            dusted_file => {
                             is => 'Text',
                             doc => 'the output fasta dusted file',
                         },
        ],
};

sub create {
    my $class = shift;
    
    my $self = $class->SUPER::create(@_);

    return $self;
}

sub execute {
    my $self = shift;
    my $fasta_file = $self->fasta_file;
    my $dusted_file = $self->dusted_file;
    my $cmd = "dust $fasta_file > $dusted_file";
    my $rv = system("dust $fasta_file > $dusted_file");

    return 1;
}
