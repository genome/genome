package Genome::Model::Tools::SmrtAnalysis::ReferenceCcs;

use strict;
use warnings;

use Genome;

class Genome::Model::Tools::SmrtAnalysis::ReferenceCcs {
    is => ['Genome::Model::Tools::SmrtAnalysis::Base'],
    has_input => [
        output_directory => {
            is => 'Text',
            doc => 'output directory',
        },
        reference_directory => {
            is => 'Text',
            doc => 'path to refRepositoryPath',
        },
        aligned_cmp_hdf5_file => {
            is => 'Text',
            doc => 'intial alignments file in cmp.h5 format',
        },
        raw_read_fasta_file => {
            is => 'Text',
            doc => 'raw read fasta file',
        },
    ],
    has_optional_input => [
        use_aln_ref_seq => {
            is => 'Boolean',
            doc => 'using the reference stored in the cmp.h5 file',
            default_value => 1,
        },
        prefix => {
            is => 'Text',
            doc => 'prefix of the output fasta and fastq files',
            default_value => 'RCCS',
        },
        min_subread_accuracy => {
            is => 'Number',
            doc => 'minimum average alignment score for sub-read identificaiton',
            default_value => 0.75,
        },
        min_fraction => {
            is => 'Number',
            doc => 'minimum fraction of the alignment length to the reference sequence length to be considered for consensus calling',
            default_value => 0,
        },
        output_cmp_hdf5_file => {
            is => 'Text',
            doc => 'output cmp.h5 file',
        },
        reference_sequence => {
            is => 'Text',
            doc => 'ref sequence fasta file',
        },
    ],
    has_optional_output => [
        frequency_count_gff_file => { },
        consensus_fastq_file => { },
        consensus_fasta_file => { },
        rccs_per_base_info_file => { },
        rccs_bam_file => { },
    ],
};

sub execute {
    my $self = shift;
    my $cmd = $self->analysis_bin .'/referenceCCSCall.py';
    if ($self->use_aln_ref_seq) {
        $cmd .= ' --useAlnRefSeq';
    }
    if ($self->min_subread_accuracy) {
        $cmd .= ' --threshold='. $self->min_subread_accuracy;
    }
    if ($self->min_fraction) {
        $cmd .= ' --minFrac='. $self->min_fraction;
    }
    my @input_files;
    my @output_files;
    if (defined($self->raw_read_fasta_file)) {
        $cmd .= ' --rawReadFasta='. $self->raw_read_fasta_file;
        push @input_files, $self->raw_read_fasta_file;
    }
    if (defined($self->reference_directory)) {
        $cmd .= ' --refRepository='. $self->reference_directory;
    }
    if (defined($self->aligned_cmp_hdf5_file)) {
        $cmd .= ' --alnCMPH5='. $self->aligned_cmp_hdf5_file;
        push @input_files, $self->aligned_cmp_hdf5_file;
    }
    if (defined($self->output_directory)) {
        $cmd .= ' --outdir='. $self->output_directory;
    }
    if (defined($self->prefix)) {
        $cmd .= ' --prefix='. $self->prefix;
    }
    if (defined($self->output_cmp_hdf5_file)) {
        $cmd .= ' --outputCMPH5='. $self->output_cmp_hdf5_file;
        push @output_files, $self->output_cmp_hdf5_file;
    }
    if (defined($self->reference_sequence)) {
        $cmd .= ' --refSeq='. $self->reference_sequence;
        push @input_files, $self->reference_sequence;
    }
    unless (defined($self->rccs_bam_file)) {
        $self->rccs_bam_file($self->output_directory .'/'. $self->prefix .'.sam');
        push @output_files, $self->rccs_bam_file;
    } else {
        die('Do not define rccs_bam_file!');
    }
    unless (defined($self->rccs_per_base_info_file)) {
        $self->rccs_per_base_info_file($self->output_directory .'/'. $self->prefix .'_pi.gz');
        push @output_files, $self->rccs_per_base_info_file;
    } else {
        die('Do not define rccs_per_base_info_file!');
    }
    unless (defined($self->consensus_fasta_file)) {
        $self->consensus_fasta_file($self->output_directory .'/'. $self->prefix .'.fasta');
        push @output_files, $self->consensus_fasta_file;
    } else {
        die('Do not define consensus_fasta_file!');
    }
    unless (defined($self->consensus_fastq_file)) {
        $self->consensus_fastq_file($self->output_directory .'/'. $self->prefix .'.fastq');
        push @output_files, $self->consensus_fastq_file;
    } else {
        die('Do not define consensus_fastq_file!');
    }
    unless (defined($self->frequency_count_gff_file)) {
        $self->frequency_count_gff_file($self->output_directory .'/'. $self->prefix .'_count.gff');
        push @output_files, $self->frequency_count_gff_file;
    } else {
        die('Do not define frequency_count_gff_file!');
    }
    $self->shellcmd(
        cmd => $cmd,
        input_files => \@input_files,
        output_files => \@output_files,
        skip_if_output_is_present => 0,
    );
    return 1;
}


1;
