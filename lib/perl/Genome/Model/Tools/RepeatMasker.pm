package Genome::Model::Tools::RepeatMasker;

use strict;
use warnings;

use Genome;
use File::Basename;

class Genome::Model::Tools::RepeatMasker {
    is => ['Command'],
    is_abstract => 1,
    has_input => [
        fasta_file => {
            is => 'Text',
            doc => 'The fasta file of sequences to mask',
        },
        output_directory => {
            is => 'Text',
            doc => 'The output directory where files are written',
        },
        sensitivity => {
            is => 'Text',
            default_value => '',
            doc => 'the sensitivity level to use which effects speed',
            valid_values => ['-s','-q','-qq',''],
        },
        species => {
            is => 'Text',
            doc => 'the species specific library to use',
            default_value => 'human',
        },
        mask => {
            is => 'Text',
            doc => 'the type of masking to perform on the masked file, RepeatMasker default is Ns',
            default_value => '-n',
            valid_values => ['-small','-xsmall','-x','-n'],
        },
    ],
};

1;
