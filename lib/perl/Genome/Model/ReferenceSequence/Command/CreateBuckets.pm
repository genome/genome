package Genome::Model::ReferenceSequence::Command::CreateBuckets;

use strict;
use warnings;

use Genome;

class Genome::Model::ReferenceSequence::Command::CreateBuckets {
    is => 'Genome::Command::DelegatesToResult',
    has => [
        reference_sequence_build => {
            is => 'Genome::Model::Build::ReferenceSequence',
            doc => 'The reference to bucketize',
        },
    ],
};

sub result_class { return 'Genome::Model::Build::ReferenceSequence::Buckets'; }

1;
