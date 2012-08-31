package Genome::Model::Tools::Fasta::To::Scf;

use strict;
use warnings;

use Genome;

use Bio::Seq::SequenceTrace;

class Genome::Model::Tools::Fasta::To::Scf {
    is           => 'Genome::Model::Tools::Fasta::To',
    has_optional => [
        dir  => {
            is      => 'String',
            doc     => 'The path of chromat_dir to place scf trace files',
            default => './chromat_dir',
        },
    ],
};


sub help_brief {
    "convert fasta sequence (and qual files) to scf trace files"
}


sub help_detail {                           
    return <<EOS 
This tool convert fasta sequence (and quality values) to scf trace files using BioPerl.
All generated scf files are named after their corrseponding fasta id without any suffix.
EOS
}


sub _format_type {
    return 'scf';
}


sub _param_type {
    return '-target';
}


sub _param_value {
    return Bio::Seq::SequenceTrace->new(-swq => pop);
}





1;
