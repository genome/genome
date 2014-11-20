package Genome::VariantReporting::Suite::Joinx::Homopolymer::Run;

use strict;
use warnings FATAL => 'all';
use Genome;

class Genome::VariantReporting::Suite::Joinx::Homopolymer::Run {
    is => 'Genome::VariantReporting::Framework::Component::Expert::Command',
    has_input => [
        joinx_version => {
            is  => 'Version',
            doc => 'joinx version to use',
        },
        homopolymer_list_id => {
            is  => 'Text',
            doc => 'Homopolymer bed file feature list id',
        },
        max_length => {
            is  => 'Integer',
            doc => 'maximum indel length to annotate as in the homopolymer, default is 2',
        },
        info_string => {
            is  => 'Text',
            doc => 'name of per-allele info field to store the annotation, default is HOMP_FILTER',
        }
    ],
};

sub name {
    'homopolymer';
}

sub result_class {
    'Genome::VariantReporting::Suite::Joinx::Homopolymer::RunResult';
}
