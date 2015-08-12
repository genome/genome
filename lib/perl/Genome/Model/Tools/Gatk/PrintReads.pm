package Genome::Model::Tools::Gatk::PrintReads;

use strict;
use warnings;

use above 'Genome';

require File::Basename;

use Genome::Model::Tools::Gatk::WithNumberOfCpuThreads;

class Genome::Model::Tools::Gatk::PrintReads {
    doc => "Run GATK with the 'PrintReads' tool",
    is => [qw/ Genome::Model::Tools::Gatk::Base Genome::Model::Tools::Gatk::WithNumberOfCpuThreads /],
    has => {
        input_bams => {
            is_input => 1,
            is_many => 1,
            is => 'Text',
            gatk_param_name => '-I',
            doc => 'The path to the original bam(s) you would like to assess',
        },
        reference_fasta => {
            is_input => 1,
            is => 'Text',
            gatk_param_name => '-R',
            doc => "The path to the reference fasta you would like to run against" ,
        },
        output_bam => {
            is_input => 1,
            is_output => 1,
            is => 'Text',
            gatk_param_name => '-o',
            doc => 'The path to where you would like to create the output bam file',
        },
        bqsr => {
            is_input => 1,
            is => 'Text',
            is_optional => 1,
            gatk_param_name => '--BQSR',
            doc => 'The input covariates table file which enables on-the-fly base quality score recalibration',
        },
        downsample_coverage => {
            is_input => 1,
            is_optional => 1,
            is => 'Number',
            gatk_param_name => '--downsample_coverage',
            doc => 'Downsample BAM to desired coverage',
        },
        number => {
            is_input => 1,
            is_optional => 1,
            is => 'Integer',
            gatk_param_name => '--number',
            doc => 'Print the first n reads from the file, discarding the rest',
        },
        platform => {
            is_input => 1,
            is_optional => 1,
            is => 'Text',
            gatk_param_name => '--platform',
            doc => 'Exclude all reads with this platform from the output',
        },
        read_group => {
            is_input => 1,
            is_optional => 1,
            is => 'Text',
            gatk_param_name => '--readGroup',
            doc => 'Exclude all reads with this read group from the output',
        },
        sample_file => {
            is_input => 1,
            is_many => 1,
            is_optional => 1,
            is => 'Text',
            gatk_param_name => '--sample_file',
            doc => 'File containing a list of samples (one per line). Can be specified multiple times',
        },
        sample_name => {
            is_input => 1,
            is_many => 1,
            is_optional => 1,
            is => 'Text',
            gatk_param_name => '--sample_name',
            doc => 'Sample name to be included in the analysis. Can be specified multiple times.',
        },
        simplify => {
            is_input => 1,
            is_optional => 1,
            is => 'Boolean',
            gatk_param_name => '--simplify',
            doc => 'Simplify all reads. Erase all extra attributes in the read but keep the read group information.',
        },
    },
};

sub help_brief {
    "Run GATK with the 'PrintReads' tool"
}

sub help_synopsis {
    return <<EOS
    gmt gatk print-reads --input-bams input1.bam,input2.bam --reference-fasta ref.fasta --output-bam output.bam --bqsr recalibration_report.grp
EOS
}

sub analysis_type {
    return 'PrintReads';
}

sub _shellcmd_extra_params {
    my $self = shift;

    my @input_files = $self->input_bams, $self->reference_fasta, $self->sample_file;
    push @input_files, $self->bqsr if $self->bqsr;

    return (
        input_files => \@input_files,
        output_files => [$self->output_bam],
    );
}

sub execute {
    my $self = shift;

    my $rv = $self->SUPER::_execute_body;

    # Rename the bam index. It gets named without a .bam in the name
    my ($bam_basename, $bam_dirname) = File::Basename::fileparse($self->output_bam);
    $bam_basename =~ s/\.bam$//;
    Genome::Sys->rename($bam_dirname.'/'.$bam_basename.'.bai', $bam_dirname.'/'.$bam_basename.'.bam.bai');

    return $rv;
}

1;
