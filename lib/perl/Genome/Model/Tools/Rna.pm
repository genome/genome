package Genome::Model::Tools::Rna;

use strict;
use warnings;

use Genome;

class Genome::Model::Tools::Rna {
    is => ['Command'],
};

sub help_brief {
    "Tool to plot output from Tophat and Cufflinks"
}

sub help_detail {
    "Takes junctions.bed and transcript.gtf files to plot assembled transcript
        and splice junctions for a specified genomic location"
}

1;
