package Genome::Model::Tools::Varscan::Consensus;

use strict;
use warnings;

use above 'Genome';

class Genome::Model::Tools::Varscan::Consensus {
    is => 'Genome::Model::Tools::Varscan',

    has_input => [
        bam_file => {
            is => 'Path',
            doc => "Path to BAM file",
        },
        output_file => {
            is => 'Path',
            is_output => 1,
            doc => "Path to output file",
        },
        ref_fasta => {
            is => 'Path',
            example_values => ['/gscmnt/sata420/info/model_data/2857786885/build102671028/all_sequences.fa'],
            doc => "Reference FASTA file for BAMs",
        },
    ],
    has_optional_input => [
        min_coverage => {
            is => 'Number',
            default => 3,
            doc => "Minimum base coverage to report readcounts",
        },
        min_avg_qual => {
            is => 'Number',
            default => 20,
            doc => "Minimum base quality to count a read",
        },
        min_var_freq => {
            is => 'Number',
            default => 0.20,
            doc => "Minimum variant allele frequency to call a variant",
        },
        vcf_sample_name => {
            is => 'Text',
            doc => 'If set, and --output-vcf is set, this name will be used instead of Sample1',
        },
        output_vcf => {
            is => 'Boolean',
            default => 0,
            doc => "If set to 1, tells VarScan to output in VCF format (rather than native CNS)",
        },
        position_list_file => {
            is => 'Path',
            doc => "Optionally, provide a tab-delimited list of positions to be given to SAMtools with -l",
        },
    ],
};

sub sub_command_sort_position { 12 }

sub help_brief {
    "Run VarScan consensus calling for one BAM file"
}

sub help_synopsis {
    return <<EOS
Runs mpileup and then VarScan consensus calling (pileup2cns) on a single BAM file
EXAMPLE:    gmt varscan consensus --bam-file sample.bam --ref-fasta reference.fa --output sample.bam.varScan.cns
EOS
}

sub help_detail {
    return <<EOS

EOS
}

sub execute {
    my $self = shift;

    Genome::Sys->shellcmd(
        cmd => $self->cmd,
        input_files => [
            $self->input_files,
        ],
        output_files => [
            $self->output_file,
        ],
    );

    return 1;
}

sub cmd {
    my $self = shift;

    my @cmd = (
        $self->varscan_command, sprintf('<(%s)', $self->samtools_command),
        '--min-coverage', $self->min_coverage,
        '--min-var-freq', $self->min_var_freq,
        '--min-avg-qual', $self->min_avg_qual,
        $self->output_vcf_string,
        sprintf("> %s 2>/dev/null", $self->output_file), # is it wise to throw away stderr?
    );
    return $self->java_command_line(join(' ', @cmd));
}

sub varscan_command {
    my $self = shift;
    if ($self->output_vcf) {
        return 'mpileup2cns';
    } else {
        return 'pileup2cns';
    }
}

sub samtools_command {
    my $self = shift;

    my @cmd = (
        $self->samtools_path,
        'mpileup', '-B',
        '-f', $self->ref_fasta,
        '-q', '10',
        $self->position_list_file_string,
        $self->bam_file,
    );

    return join(' ', @cmd);
}

sub output_vcf_string {
    my $self = shift;
    if ($self->output_vcf) {
        my $str = '--output-vcf 1';
        if ($self->vcf_sample_name) {
            $str .= sprintf(' --vcf-sample-list <(echo "%s")', $self->vcf_sample_name);
        }
        return $str;
    } else {
        return '';
    }
}

sub position_list_file_string {
    my $self = shift;

    if($self->position_list_file) {
        return '-l ' . $self->position_list_file;
    } else {
        return '';
    }
}

sub input_files {
    my $self = shift;

    my @files = ($self->ref_fasta, $self->bam_file);
    if ($self->position_list_file) {
        push @files, $self->position_list_file;
    }
    return @files;
}


1;
