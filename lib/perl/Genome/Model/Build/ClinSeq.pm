package Genome::Model::Build::ClinSeq;
use strict;
use warnings;

use Genome;
class Genome::Model::Build::ClinSeq {
  is => ['Genome::Model::Build',
    'Genome::Model::Build::ClinSeq::FileAccessors',
    'Genome::Model::Build::ClinSeq::InputBuilds'],
};

1;
