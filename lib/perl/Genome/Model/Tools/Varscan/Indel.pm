package Genome::Model::Tools::Varscan::Indel;

use strict;
use warnings;

use above 'Genome';

class Genome::Model::Tools::Varscan::Indel {
    is => 'Genome::Model::Tools::Varscan::Basic',
};

sub help_brief {
    "Run VarScan indel calling for one BAM file"
}

sub help_synopsis {
    return <<EOS
Runs mpileup and then VarScan indel calling (pileup2indel) on a single BAM file
EXAMPLE:    gmt varscan indel --bam-file sample.bam --ref-fasta reference.fa --output sample.bam.varScan.cns
EOS
}

sub help_detail {
    return <<EOS

EOS
}

sub varscan_command {
    my $self = shift;

    # http://www.biostars.org/p/6259/ says not to use pileup for indels,
    # so we have it set to use mpileup.
    if ($self->output_vcf) {
        return 'mpileup2indel';
    } else {
        return 'mpileup2indel';
    }
}

1;
