package Genome::Model::Tools::Picard::CreateSequenceDictionary;

use strict;
use warnings;

use Genome;
use File::Basename;


class Genome::Model::Tools::Picard::CreateSequenceDictionary {
    is => 'Genome::Model::Tools::Picard::Base',

    has_input => [
        reference_fasta => {
            doc => 'Input reference fasta or fasta.gz',
            picard_param_name => 'REFERENCE',
        },
        species => {
            is_optional => 1,
            doc => 'Put into SP field of sequence dictionary entry.',
            picard_param_name => 'SPECIES',
        },
        genome_assembly => {
            is_optional => 1,
            doc => 'Put into AS field of sequence dictionary entry if supplied',
            picard_param_name => 'GENOME_ASSEMBLY',
        },
        uri => {
            is_optional => 1,
            doc => 'Put into UR field of sequence dictionary entry. If not supplied, input reference file is used.',
            picard_param_name => 'URI',
        },
        output_file => {
            is => 'String',
            doc => 'Output SAM or BAM file containing only the sequence dictionary.',
            picard_param_name => 'OUTPUT',
        },
        truncate_names_at_whitespace => {
            is => 'Boolean',
            default_value => 1,
            is_optional => 1,
            doc => 'Make sequence name the first word from the > line in the fasta file. By default the entire contents of the > line is used, excluding leading and trailing whitespace.',
            picard_param_name => 'TRUNCATE_NAMES_AT_WHITESPACE',
        },
        num_sequences => {
            is => 'Number',
            is_optional => 1,
            doc => '	Stop after writing this many sequences. For testing. Default value: 2147483647.',
            picard_param_name => 'NUM_SEQUENCES',
        },
    ],
};

sub help_brief {
    'Read fasta or fasta.gz containing reference sequences, and write as a SAM or BAM file with only sequence dictionary.';
}

sub help_detail {
    return <<EOS
    For Picard documentation of this command see:
    http://broadinstitute.github.io/picard/command-line-overview.html#CreateSequenceDictionary
EOS
}

sub _jar_name {
    return 'CreateSequenceDictionary.jar';
}

sub _java_class {
    return qw(picard sam CreateSequenceDictionary);
}

1;
