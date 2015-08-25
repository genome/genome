package Genome::Model::Tools::Gatk::RealignerTargetCreator;

use strict;
use warnings;
use Genome;

class Genome::Model::Tools::Gatk::RealignerTargetCreator {
    doc => "Run GATK with the 'RealignerTargetCreator' tool",
    is => [qw/ Genome::Model::Tools::Gatk::Base Genome::Model::Tools::Gatk::WithNumberOfThreads /],
    has_input => [
        known => {
            is => 'Text',
            doc => 'The file of known indels',
            is_optional => 1,
            is_many => 1,
            gatk_param_name => '--known',
        },
        input_bam => {
            is => 'Text',
            doc => 'The path to the original bam you would like to be realigned',
            gatk_param_name => '-I',
        },
        reference_fasta => {
            is => 'Text',
            doc => "Reference Fasta" ,
            gatk_param_name => '-R',
        },
        output_intervals => {
            is => 'Text',
            doc => "File of intervals to target for realignment",
            gatk_param_name => '-o',
        },
    ],
};

sub help_brief {
    "Run GATK with the 'RealignerTargetCreator' tool"
}

sub help_synopsis {
    return <<EOS
    gmt gatk realigner-indels-target-creator --known indels.vcf --input-bam my_existing.bam --reference-fasta my.fa --output-intervals out.intervals
EOS
}

sub analysis_type {
    return 'RealignerTargetCreator';
}

sub _shellcmd_extra_params {
    my $self = shift;

    return (
       input_files => [$self->input_bam, $self->reference_fasta, $self->known],
       output_files => [$self->output_intervals],
    );
}

1;
