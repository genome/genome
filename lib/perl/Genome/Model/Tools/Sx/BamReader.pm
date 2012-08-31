package Genome::Model::Tools::Sx::BamReader;

use strict;
use warnings;

use Genome;

require File::Which;

class Genome::Model::Tools::Sx::BamReader {
    is => 'Genome::Model::Tools::Sx::SamReader',
};

sub _cmd_for_file {
    my ($self, $file) = @_;

    my $executable = File::Which::which('samtools');
    if ( not $executable ) {
        die "Failed to find executable 'samtools'. Cannot read from BAM file: $file";
    }

    return'samtools view '.$file.' |';
}

1;

