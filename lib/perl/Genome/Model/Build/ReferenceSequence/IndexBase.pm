package Genome::Model::Build::ReferenceSequence::IndexBase;

use strict;
use warnings;

use Genome;

class Genome::Model::Build::ReferenceSequence::IndexBase {
    is_abstract => 1,
    is => ['Genome::SoftwareResult::Stageable', 'Genome::SoftwareResult::WithNestedResults'],
};

1;
