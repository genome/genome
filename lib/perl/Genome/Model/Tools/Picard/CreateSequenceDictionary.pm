package Genome::Model::Tools::Picard::CreateSequenceDictionary;

use strict;
use warnings;

use Genome;
use File::Basename;


class Genome::Model::Tools::Picard::CreateSequenceDictionary {
    is  => 'Genome::Model::Tools::Picard',
    has_input => [
        reference_fasta => {
            doc => 'Input reference fasta or fasta.gz',
        },
        species => {
            is_optional => 1,
            doc => 'Put into SP field of sequence dictionary entry.',
        },
        genome_assembly  => {
            is_optional => 1,
            doc => 'Put into AS field of sequence dictionary entry if supplied',
        },
        uri => {
            is_optional => 1,
            doc => 'Put into UR field of sequence dictionary entry. If not supplied, input reference file is used.',
        },
        output_file => {
            is  => 'String',
            doc => 'Output SAM or BAM file containing only the sequence dictionary.',
        },
        truncate_names_at_whitespace => {
            is => 'Boolean',
            default_value => 1,
            is_optional => 1,
            doc => 'Make sequence name the first word from the > line in the fasta file. By default the entire contents of the > line is used, excluding leading and trailing whitespace.',
        },
        num_sequences => {
            is => 'Number',
            is_optional => 1,
            doc => '	Stop after writing this many sequences. For testing. Default value: 2147483647.',
        },
    ],
};

sub help_brief {
    'Read fasta or fasta.gz containing reference sequences, and write as a SAM or BAM file with only sequence dictionary.';
}

sub help_detail {
    return <<EOS
    For Picard documentation of this command see:
    http://picard.sourceforge.net/command-line-overview.shtml#CreateSequenceDictionary
EOS
}

sub execute {
    my $self = shift;

    my $cmd = $self->picard_path .'/CreateSequenceDictionary.jar net.sf.picard.sam.CreateSequenceDictionary';
    $cmd   .= ' OUTPUT='. $self->output_file  .' REFERENCE='. $self->reference_fasta;
    if (defined($self->species)) {
        $cmd .= ' SPECIES='. $self->species;
    }
    if (defined($self->genome_assembly)) {
        $cmd .= ' GENOME_ASSEMBLY='. $self->genome_assembly;
    }
    if (defined($self->uri)) {
        $cmd .= ' URI='. $self->uri;
    }
    if (defined($self->num_sequences)) {
        $cmd .= ' NUM_SEQUENCES='. $self->num_sequences;
    }
    unless ($self->truncate_names_at_whitespace) {
        $cmd .= ' TRUNCATE_NAMES_AT_WHITESPACE=false';
    }
    $self->run_java_vm(
        cmd          => $cmd,
        input_files  => [$self->reference_fasta],
        output_files => [$self->output_file],
        skip_if_output_is_present => 0,
    );
    return 1;
}


1;
