package Genome::Model::Tools::Fasta::To::Phd;

use strict;
use warnings;

use Genome;

class Genome::Model::Tools::Fasta::To::Phd {
    is           => 'Genome::Model::Tools::Fasta::To',
    has_optional => [
        dir  => {
            is      => 'String',
            doc     => 'The path of phd_dir to place phd files',
            default => './phd_dir',
        },
        time => {
            is      => 'String',
            doc     => 'time stamp inside phd file, often need sync with timestamp in acefile',
        }
    ],
};


sub help_brief {
    "convert fasta sequence (and qual files) to phd files"
}


sub help_detail {                           
    return <<EOS 
This tool convert fasta sequence (and quality values) to phd files using BioPerl. All 
generated phd files get '.phd.1' as suffix. Trace index (the third column) is made up
by increment 10 and starting from 0.
EOS
}


sub _format_type {
    return 'phd';
}

1;
