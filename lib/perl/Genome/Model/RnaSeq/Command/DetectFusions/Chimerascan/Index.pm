package Genome::Model::RnaSeq::Command::DetectFusions::Chimerascan::Index;

use strict;
use warnings;

use Genome;

class Genome::Model::RnaSeq::Command::DetectFusions::Chimerascan::Index {
    is => 'Genome::Model::RnaSeq::Command::DetectFusions::Chimerascan::IndexBase',
    doc => 'create the annotation index used inside chimerascan',
};

sub help_synopsis {
    return <<EOS
 genome model rna-seq detect-fusions chimerascan-index --version=1.2.3 --bowtie-version=1.2.3 --reference-build=1234 --annotation-build=1234

EOS
}

sub help_detail {
    return <<EOS
Create the annotation index used inside chimerascan
EOS
}

sub result_class_name {
    "Genome::Model::RnaSeq::DetectFusionsResult::Chimerascan::Index";
}

1;
