package Genome::Model::ClinSeq::Command::Converge::Base;
use strict;
use warnings;
use Genome;

class Genome::Model::ClinSeq::Command::Converge::Base {
    is => 'Command::V2',
    is_abstract => 1,
    has_input => [
        builds => { is => 'Genome::Model::Build::ClinSeq',
                    doc => 'clinseq builds to converge' },
    ],
    doc => 'converge various data types across clinseq inputs'
};

1;

