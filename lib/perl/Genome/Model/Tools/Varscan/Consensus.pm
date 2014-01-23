package Genome::Model::Tools::Varscan::Consensus;

use strict;
use warnings;

use above 'Genome';

class Genome::Model::Tools::Varscan::Consensus {
    is => 'Genome::Model::Tools::Varscan::Basic',
};

sub help_brief {
    "Run VarScan consensus calling for one BAM file"
}

sub help_synopsis {
    return <<EOS
Runs mpileup and then VarScan consensus calling (pileup2cns) on a single BAM file
EXAMPLE:    gmt varscan consensus --bam-file sample.bam --ref-fasta reference.fa --output sample.bam.varScan.cns
EOS
}

sub help_detail {
    return <<EOS

EOS
}

sub varscan_command {
    my $self = shift;
    if ($self->output_vcf) {
        return 'mpileup2cns';
    } else {
        return 'pileup2cns';
    }
}

1;
