package Genome::Model::RnaSeq::Command::DetectFusions::Chimerascan::Main;

use strict;
use warnings;

use Genome;

class Genome::Model::RnaSeq::Command::DetectFusions::Chimerascan::Main {
    is => 'Genome::Model::RnaSeq::Command::DetectFusions::Chimerascan::MainBase',
    doc => 'run the chimerascan transcript fusion detector',
};

sub help_synopsis {
    return <<EOS
 genome model rna-seq detect-fusions chimerascan index_dir/ f1.fastq f2.fast2 output_dir/

 genome model rna-seq detect-fusions chimerascan index_dir/ f1.fastq f2.fast2 output_dir/ --use-version 1.2.3 --params "-a -b -c"
EOS
}

sub help_detail {
    return <<EOS
Run the chimerascan gene fusion detector.

It is used by the RNASeq pipeline to perform fusion detection when the fusion detection strategy is set to something like:
 'chimerascan 1.2.3'

EOS
}

sub result_class {
    return "Genome::Model::RnaSeq::DetectFusionsResult::Chimerascan::VariableReadLength::Result";
}

1;
