package Genome::Model::RnaSeq::Command::DetectFusions::Chimerascan::VariableReadLength::Index;

use strict;
use warnings;

use Genome;

class Genome::Model::RnaSeq::Command::DetectFusions::Chimerascan::VariableReadLength::Index {
    is => 'Genome::Model::RnaSeq::Command::DetectFusions::Chimerascan::IndexBase',
    doc => 'create the annotation index used inside chimerascan-vrl',
};

sub help_synopsis {
    return <<EOS
 genome model rna-seq detect-fusions chimerascan-vrl-index --version=0.4.6 --bowtie-version=0.12.7 --picard-version 1.2.3 --reference-build=1234 --annotation-build=1234

EOS
}

sub help_detail {
    return <<EOS
Create the annotation index used inside chimerascan-vrl
EOS
}

sub result_class_name {
    "Genome::Model::RnaSeq::DetectFusionsResult::Chimerascan::VariableReadLength::Index";
}

1;
