package Genome::Model::RnaSeq::Command::DetectFusions::Chimerascan::VariableReadLength::Detector;

use strict;
use warnings;

use Genome;

class Genome::Model::RnaSeq::Command::DetectFusions::Chimerascan::VariableReadLength::Detector {
    is => 'Genome::Model::RnaSeq::Command::DetectFusions::Chimerascan::DetectorBase',
    doc => 'run the chimerascan-vrl (variable-read-length) transcript fusion detector',
};

sub help_synopsis {
    return <<EOS
 genome model rna-seq detect-fusions chimerascan-vrl build-id=1234
EOS
}

sub help_detail {
    return <<EOS
Run the chimerascan-vrl gene fusion detector. ChimerascanVrl extends Chimerascan to include support for
reads with varying numbers of reads.

It is used by the RNASeq pipeline to perform fusion detection when the fusion detection strategy is set to something like:
 'chimerascan-vrl 0.4.6 [-p 2 --reuse-bam=1 --bowtie-version=0.12.7]'

EOS
}

sub result_class {
    return "Genome::Model::RnaSeq::DetectFusionsResult::Chimerascan::VariableReadLength::Result";
}

1;
