package Genome::Model::Tools::Picard::SortVcf;

use strict;
use warnings;

use Genome;

my $DEFAULT_PICARD_VERSION = '1.123';

class Genome::Model::Tools::Picard::SortVcf {
    is  => 'Genome::Model::Tools::Picard::Base',
    has_input => [
        input_vcf => {
            picard_param_name => 'INPUT',
            is  => 'String',
            is_many => 1,
            doc => 'Input VCF(s) to be sorted. Multiple inputs must have the same sample names (in order) This option may be specified 0 or more times.',
        },
        output_vcf => {
            picard_param_name => 'OUTPUT',
            is  => 'String',
            doc => 'Output VCF to be written.',
        },
        sequence_dictionary => {
            picard_param_name => 'SEQUENCE_DICTIONARY',
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

sub minimum_version_required { '1.122'; }
sub _jar_name { 'SortVcf.jar'; }
sub _java_class { 'picard.vcf.SortVcf'; }

sub _shellcmd_extra_params {
    my $self = shift;
    return (
        input_files => [ $self->input_vcf ],
        output_files => [ $self->output_vcf ],
        skip_if_output_is_present => 0,
    );
}

1;

