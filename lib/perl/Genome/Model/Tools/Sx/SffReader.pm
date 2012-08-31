package Genome::Model::Tools::Sx::SffReader;

use strict;
use warnings;

use Genome;

require File::Which;

class Genome::Model::Tools::Sx::SffReader {
    is => 'Genome::Model::Tools::Sx::FastqReader',
};

sub _cmd_for_file {
    my ($self, $file) = @_;

    my $executable = File::Which::which('sff2fastq');
    if ( not $executable ) {
        die "Failed to find executable 'sff2fastq'.\nCannot read from SFF file: $file";
    }

    return "sff2fastq $file |";
}

1;

