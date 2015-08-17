package Genome::Model::Tools::Gatk::BaseRecalibrator;

use strict;
use warnings;
use Genome;

class Genome::Model::Tools::Gatk::BaseRecalibrator {
    doc => "Run GATK with the 'BaseRecalibrator' tool",
    is => [qw/ Genome::Model::Tools::Gatk::Base Genome::Model::Tools::Gatk::WithNumberOfCpuThreads /],
    has => [
        input_bam => {
            is_input => 1,
            is => 'Text',
            gatk_param_name => '-I',
            doc => 'The path to the original bam you would like to assess',
        },
        reference_fasta => {
            is_input => 1,
            is => 'Text',
            gatk_param_name => '-R',
            doc => "The path to the reference fasta you would like to run against" ,
        },
        known_sites => {
            is_input => 1,
            is_many => 1,
            is => 'Text',
            gatk_param_name => '-knownSites',
            doc => 'A database of known polymorphic sites to skip over in the recalibration algorithm',
        },
        output_recalibration_table => {
            is_input => 1,
            is_output => 1,
            is => 'Text',
            gatk_param_name => '-o',
            doc => 'The path to where you would like to create the output recalibration table file',
        },
    ],
};

sub help_brief {
    "Run GATK with the 'BaseRecalibrator' tool"
}

sub help_synopsis {
    return <<EOS
    gmt gatk base-recalibator --input-bam my_reads.bam --reference-fasta my.fa --known-sites path/to/sites.vcf,path/to/more/sites.vcf --output-recalibration-table recal_data.grp
EOS
}

sub analysis_type {
    return 'BaseRecalibrator';
}

sub _shellcmd_extra_params {
    my $self = shift;

    return (
        input_files => [$self->input_bam, $self->reference_fasta, $self->known_sites],
        output_files => [$self->output_recalibration_table],
    );
}

1;
