package Genome::Model::Tools::Gatk::HaplotypeCaller;

use strict;
use warnings;
use Genome;

class Genome::Model::Tools::Gatk::HaplotypeCaller {
    doc => "Run GATK with the 'HaplotypeCaller' tool",
    is => 'Genome::Model::Tools::Gatk::Base',
    has => [
        input_bam => {
            is_input => 1,
            is => 'Text',
            is_many => 1,
            gatk_param_name => '-I',
            doc => 'The path to the original bam(s) you would like to analyze',
        },
        reference_fasta => {
            is_input => 1,
            is => 'Text',
            gatk_param_name => '-R',
            doc => "The path to the reference fasta you would like to run against" ,
        },
        dbsnp_vcf => {
            is_input => 1,
            is => 'Text',
            is_optional => 1,
            gatk_param_name => '--dbsnp',
            doc => 'The dbSNP file to use to annotate the output',
        },
        output_vcf => {
            is_input => 1,
            is_output => 1,
            is => 'Text',
            gatk_param_name => '-o',
            doc => 'The path to where you would like to create the output VCF file',
        },
        emit_reference_confidence => {
            is_input => 1,
            is => 'Text',
            default_value => 'NONE',
            valid_values => ['NONE', 'GVCF', 'BP_RESOLUTION'],
            gatk_param_name => '-ERC',
            doc => 'Whether the trimming intervals will be used to emit reference confidence',
        },
        gvcf_gq_bands => {
            is_input => 1,
            is => 'Number',
            is_many => 1,
            gatk_param_name => '-GQB',
            doc => 'GQ thresholds for reference confidence bands',
            is_optional => 1,
        },
        intervals => {
            is_input => 1,
            is => 'Text',
            gatk_param_name => '-L',
            doc => 'restrict run to these regions (either a file or explicit intervals)',
            is_many => 1,
            is_optional => 1,
        },
        read_filters => {
            is => 'Text',
            is_many => 1,
            is_optional => 1,
            doc => 'Filters to apply to reads before analysis.',
            gatk_param_name => '-rf',
        },
    ],
};

sub help_brief {
    "Run GATK with the 'HaplotypeCaller' tool"
}

sub help_synopsis {
    return <<EOS
    gmt gatk haplotype-caller --input-bam my_reads.bam --reference-fasta my.fa --dbsnp-vcf path/to/dbsnp.vcf --output-vcf output.g.vcf --emit-reference-confidence GVCF
EOS
}

sub analysis_type {
    return 'HaplotypeCaller';
}

sub _shellcmd_extra_params {
    my $self = shift;

    my @inputs = $self->input_bam, $self->reference_fasta;
    push @inputs, $self->dbsnp_vcf if $self->dbsnp_vcf;

    return (
        input_files => \@inputs,
        output_files => [$self->output_vcf],
    );
}

sub _cmdline_args {
    my $self = shift;

    my @args = $self->SUPER::_cmdline_args(@_);

    if ($self->version eq '3.4' and
            $self->output_vcf =~ /g\.vcf\.gz$/ and
            $self->emit_reference_confidence eq 'GVCF') {
        #Handle a bug in GATK 3.4 by including these deprecated parameters.
        push @args, qw(-variant_index_type LINEAR -variant_index_parameter 128000);
    }

    return @args;
}

1;
