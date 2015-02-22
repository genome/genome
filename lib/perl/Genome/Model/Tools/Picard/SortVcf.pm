package Genome::Model::Tools::Picard::SortVcf;

use strict;
use warnings;

use Genome;

my $DEFAULT_PICARD_VERSION = '1.123';

class Genome::Model::Tools::Picard::SortVcf {
    is  => 'Genome::Model::Tools::Picard',
    has_input => [
        input_vcf => {
            is  => 'String',
            is_many => 1,
            doc => 'Input VCF(s) to be sorted. Multiple inputs must have the same sample names (in order) This option may be specified 0 or more times.',
        },
        output_vcf => {
            is  => 'String',
            doc => 'Output VCF to be written.',
        },
        sequence_dictionary => {
            is => 'String',
            doc => 'The path to the Picard sequence dictionary for the reference genome.',
            is_optional => 1,
        },
        use_version => {
            is => 'String',
            doc => 'Version must be 1.122 or greater.',
            default_value => $DEFAULT_PICARD_VERSION,
        },
    ],
};

sub help_brief {
    'Sorts one or more VCF files according to the order of the contigs in the header/sequence dictionary and then by coordinate.';
}

sub help_detail {
    return <<EOS

Documentation: http://broadinstitute.github.io/picard/command-line-overview.html#SortVcf

Sorts one or more VCF files according to the order of the contigs in the header/sequence dictionary and then by coordinate. Can accept an external sequence dictionary. If no external dictionary is supplied, multiple inputs' headers must have the same sequence dictionaries. Multiple inputs must have the same sample names (in order) 

EOS
}

sub execute {
    my $self = shift;


    my $sort_cmd = $self->picard_path .'/SortVcf.jar picard.vcf.SortVcf OUTPUT='.$self->output_vcf;
    my @input_files;
    for my $input_vcf ($self->input_vcf) {
        $sort_cmd .= ' INPUT='. $input_vcf;
        push @input_files, $input_vcf;
    }
    if ($self->sequence_dictionary) {
        $sort_cmd .= ' SEQUENCE_DICTIONARY='. $self->sequence_dictionary;
        push @input_files, $self->sequence_dictionary;
    }
    $self->run_java_vm(
        cmd => $sort_cmd,
        input_files => \@input_files,
        output_files => [$self->output_vcf],
        skip_if_output_is_present => 0,
    );
    return 1;
}


1;
