package Genome::InstrumentData::AlignmentResult::Command::PicardRnaSeqMetrics;

use strict;
use warnings;

use Genome;
use Sys::Hostname;
use File::Path;

class Genome::InstrumentData::AlignmentResult::Command::PicardRnaSeqMetrics {
    is => ['Command::V2'],
    has_input => [
        alignment_result_id => {
            is => 'Number',
            doc => 'ID of the result for the alignment data upon which to run',
        },
        reference_build_id => {
            is => 'Number',
            doc => 'ID of the reference sequence used to generate metrics.',
        },
        annotation_build_id => {
            is => 'Number',
            doc => 'ID of the annotation used to generate metrics.',
        },
        metrics_directory => {
            is => 'Text',
            doc => 'The directory to write the Picard RnaSeq metrics.',
        },
    ],
    has_param => [
        picard_version => {
            is => 'Text',
            doc => 'The version of picard to use.',
        },
    ],
    has => [
        alignment_result => {
            is => 'Genome::SoftwareResult',
            id_by => 'alignment_result_id',
            doc => 'the alignment data upon which to run',
        },
        reference_build => {
            is => 'Genome::Model::Build::ReferenceSequence',
            id_by => 'reference_build_id',
            doc => 'the reference sequence upon which to run',
        },
        annotation_build => {
            is => 'Genome::Model::Build::ImportedAnnotation',
            id_by => 'annotation_build_id',
            doc => 'the annotation upon which to run',
        },
    ],
};

sub ribo_intervals_filename {
    return Genome::Model::Build::RnaSeq->default_ribo_intervals_filename;
}

sub mRNA_ref_flat_filename {
    return Genome::Model::Build::RnaSeq->default_mRNA_ref_flat_filename;
}

sub pie_chart_filename {
    return Genome::Model::Build::RnaSeq->default_pie_chart_filename;
}

sub metrics_output_filename {
    return Genome::Model::Build::RnaSeq->default_metrics_output_filename;
}

sub chart_output_filename {
    return Genome::Model::Build::RnaSeq->default_chart_output_filename;
}

sub execute {
    my $self = shift;
    
    my $metrics_directory = $self->metrics_directory;
    unless (-d $metrics_directory) {
        Genome::Sys->create_directory($metrics_directory);
    }
    
    my $annotation_build = $self->annotation_build;
    my $reference_build = $self->reference_build;
    my $alignment_result = $self->alignment_result;
    my $picard_version = $self->picard_version;
    
    my $bam_file = $alignment_result->bam_file;
    
    my $reference_fasta_file = $reference_build->full_consensus_path('fa');
    unless (-s $reference_fasta_file) {
        $self->error_message("Reference FASTA File ($reference_fasta_file) is missing");
        return;
    }
    
    # Resolve output file paths
    my $ribo_intervals = $metrics_directory .'/'. $self->ribo_intervals_filename;
    my $mRNA_ref_flat = $metrics_directory .'/'. $self->mRNA_ref_flat_filename;
    my $pie_chart_file = $metrics_directory .'/'. $self->pie_chart_filename;
    my $metrics_output_file = $metrics_directory .'/'. $self->metrics_output_filename;
    my $chart_output_file = $metrics_directory .'/'. $self->chart_output_filename;
    
    # Get all MT and rRNA annotation in intervals format
    my $rRNA_MT_gtf_file = $annotation_build->rRNA_MT_file('gtf',$reference_build->id,0);
    unless(-s $rRNA_MT_gtf_file) {
        $rRNA_MT_gtf_file = $annotation_build->rRNA_MT_file('gtf',undef,0);
    }
    my $seqdict_file = $reference_build->get_sequence_dictionary('sam',$reference_build->species_name,$picard_version);
    unless (Genome::Model::Tools::Gtf::ToIntervals->execute(
        gtf_file => $rRNA_MT_gtf_file,
        seqdict_file => $seqdict_file,
        interval_file => $ribo_intervals,
    )) {
        $self->error_message('Failed to convert the rRNA_MT GTF file '. $rRNA_MT_gtf_file .' to intervals: '. $ribo_intervals);
        return;
    }

    # Get all mRNA annotation in RefFlat format
    my $mRNA_gtf_file = $annotation_build->annotation_file('gtf',$reference_build->id,0);
    unless(-s $mRNA_gtf_file) {
        $mRNA_gtf_file = $annotation_build->annotation_file('gtf',undef,0);
    }
    unless (Genome::Model::Tools::Gtf::ToRefFlat->execute(
        input_gtf_file => $mRNA_gtf_file,
        output_file => $mRNA_ref_flat,
    )) {
        $self->error_message('Failed to convert the all_sequences GTF file '. $mRNA_gtf_file .'  to RefFlat: '. $mRNA_ref_flat);
        return;
    }

    # This is wasteful, but required since BAM sort order does not match our FASTA chromosome order
    my $tmp_bam_file = Genome::Sys->create_temp_file_path();
    unless (Genome::Model::Tools::Picard::ReorderSam->execute(
        input_file => $bam_file,
        output_file => $tmp_bam_file,
        reference_file => $reference_fasta_file,
        use_version => $picard_version,
    )) {
        $self->error_message('Failed to reorder BAM file: '. $bam_file);
        return;
    }
    
    unless (Genome::Model::Tools::Picard::CollectRnaSeqMetrics->execute(
        input_file => $tmp_bam_file,
        output_file => $metrics_output_file,
        refseq_file => $reference_fasta_file,
        ribosomal_intervals_file => $ribo_intervals,
        ref_flat_file => $mRNA_ref_flat,
        use_version => $picard_version,
        chart_output => $chart_output_file,
    )) {
        $self->error_message('Failed to run Picard CollectRnaSeqMetrics for alignment result: '. $alignment_result->id);
        return;
    }
    unless (Genome::Model::Tools::Picard::PlotRnaSeqMetrics->execute(
        input_file => $metrics_output_file,
        output_file => $pie_chart_file,
        label => $alignment_result->id,
    )) {
        $self->error_message('Failed to run PlotRnaSeqMetrics for alignment result: '. $alignment_result->id);
        return;
    }
    return 1;
}


1;
