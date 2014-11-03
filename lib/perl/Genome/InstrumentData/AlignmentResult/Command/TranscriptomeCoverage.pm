package Genome::InstrumentData::AlignmentResult::Command::TranscriptomeCoverage;

use strict;
use warnings;

use Genome;
use Sys::Hostname;
use File::Path;

class Genome::InstrumentData::AlignmentResult::Command::TranscriptomeCoverage {
    is => ['Command::V2'],
    has_input => [
        alignment_result_id => {
            is => 'Number',
            doc => 'ID of the result for the alignment data upon which to run coverage',
        },
        reference_build_id => {
            is => 'Number',
            doc => 'ID of the reference sequence used to generate metrics.',
        },
        annotation_build_id => {
            is => 'Number',
            doc => 'ID of the annotation used to generate metrics.',
        },
        coverage_directory => {
            is => 'Text',
            doc => 'The directory to write the coverage metrics.',
        },
    ],
    has_param => [
        annotation_file_basenames => {
            is_many => 1,
            doc => 'A list of annotation file basenames to generate coverage of',
        },
        merge_annotation_features => {
            is => 'Text',
            doc => 'The option to generate coverage on gene-level squashed exons or each unique exon.',
            valid_values => ['yes','no','both'],
        },
        mask_reference_transcripts => {
            is_optional => 1,
            doc => 'A list of transcript types to mask coverage',
            valid_values => ['rRNA','MT','pseudogene','rRNA_MT','rRNA_MT_pseudogene'],
        },
    ],
    has => [
        alignment_result => {
            is => 'Genome::SoftwareResult',
            id_by => 'alignment_result_id',
            doc => 'the alignment data upon which to run coverage',
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

sub execute {
    my $self = shift;

    # Reference inputs
    my $reference_build = $self->reference_build;
    my $reference_fasta_file = $reference_build->full_consensus_path('fa');
    die $self->error_message("Reference FASTA File ($reference_fasta_file) is missing") unless -s $reference_fasta_file;

    # Annotation inputs
    my $annotation_build = $self->annotation_build;

    # Alignment inputs
    my $bam_file = $self->alignment_result->bam_file;

    my $coverage_directory = $self->coverage_directory;
    unless (-d $coverage_directory) {
        Genome::Sys->create_directory($coverage_directory);
    }
    for my $annotation_basename ($self->annotation_file_basenames) {
        my $annotation_file_method = $annotation_basename .'_file';
        my @output_stats_files;
        if ($self->merge_annotation_features eq 'no' || $self->merge_annotation_features eq 'both') {
            my $bed_file = $annotation_build->$annotation_file_method('bed',$reference_build->id,0);
            unless ($bed_file) {
                die('Failed to find BED annotation transcripts with type: '. $annotation_basename);
            }
            if ($self->mask_reference_transcripts) {
                $bed_file = $self->remove_reference_transcripts($bed_file,0);
                unless ($bed_file) {
                    $self->error_message('Unable to remove the ROI to be masked.');
                    return;
                }
            }
            my $stats_file = $coverage_directory .'/'. $annotation_basename .'_exon_STATS.tsv';
            my $transcript_stats_file = $coverage_directory .'/'. $annotation_basename .'_transcript_STATS.tsv';
            my @output_files = ($stats_file,$transcript_stats_file);

            # TODO: Each run of ref-cov may become it's own software result...
            # TODO: Also we should run each in parallel, this can take a long time
            my %ref_cov_params = (
                alignment_file_path => $bam_file,
                roi_file_path => $bed_file,
                reference_fasta => $reference_fasta_file,
                stats_file => $stats_file,
                merged_stats_file => $transcript_stats_file,
            );
            my $refcov_cmd = Genome::Model::Tools::RefCov::RnaSeq->execute(%ref_cov_params);
            unless ($refcov_cmd and $refcov_cmd->result) {
                die('Failed to run ref_cov with params: '. %ref_cov_params);
            }
            push @output_stats_files, $stats_file;
            push @output_stats_files, $transcript_stats_file;
        }
        if ($self->merge_annotation_features eq 'yes' || $self->merge_annotation_features eq 'both') {
            my $squashed_bed_file = $annotation_build->$annotation_file_method('bed',$reference_build->id,1);
            unless ($squashed_bed_file) {
                $self->warning_message('Failed to find squashed '. $annotation_file_method .' BED file for reference build '. $reference_build->id .' in: '. $annotation_build->data_directory);
                next;
            }

            if ($self->mask_reference_transcripts) {
                $squashed_bed_file = $self->remove_reference_transcripts($squashed_bed_file,1);
                unless ($squashed_bed_file) {
                    $self->error_message('Unable to remove the ROI to be masked!');
                    return;
                }
            }
            my $squashed_stats_file = $coverage_directory .'/'. $annotation_basename .'_squashed_by_gene_STATS.tsv';
            my $genes_stats_file = $coverage_directory .'/'. $annotation_basename .'_gene_STATS.tsv';
            my %ref_cov_params = (
                alignment_file_path => $bam_file,
                roi_file_path => $squashed_bed_file,
                reference_fasta => $reference_fasta_file,
                stats_file => $squashed_stats_file,
                merged_stats_file => $genes_stats_file,
                merge_by => 'gene',
            );
            my $refcov_cmd = Genome::Model::Tools::RefCov::RnaSeq->execute(%ref_cov_params);
            unless ($refcov_cmd and $refcov_cmd->result) {
                die('Failed to run ref_cov with params: '. %ref_cov_params);
            }
            push @output_stats_files, $squashed_stats_file;
            push @output_stats_files, $genes_stats_file;
        }
        for my $stats_output_file (@output_stats_files) {
            my ($basename,$dirname,$suffix) = File::Basename::fileparse($stats_output_file,qw/\.tsv/);
            my $summary_output_file = $dirname . $basename .'.txt';
            my $stats_cmd = Genome::Model::Tools::BioSamtools::StatsSummary->execute(
                stats_file => $stats_output_file,
                output_file => $summary_output_file,
            );
            unless ($stats_cmd and $stats_cmd->result) {
                die('Failed to generate stats sumamry for stats file: '. $stats_output_file);
            }
        }
    }
    # TODO: run once using squashed transcriptome ie. merge entire BED regardless of annotation?

    return 1;
}

sub remove_reference_transcripts {
    my $self = shift;
    my $bed_file = shift;
    my $squashed_flag = shift;

    my $annotation_build = $self->annotation_build;
    my $reference_build = $self->reference_build;

    my $annotation_file_method = $self->mask_reference_transcripts .'_file';
    my $mask_bed_file = $annotation_build->$annotation_file_method('bed',$reference_build->id,$squashed_flag);
    my $tmp_bed_file = Genome::Sys->create_temp_file_path();
    my $mask_cmd = Genome::Model::Tools::BedTools::Subtract->execute(
        input_file_a => $bed_file,
        input_file_b => $mask_bed_file,
        output_file => $tmp_bed_file,
    );
    unless ($mask_cmd and $mask_cmd->result) {
        $self->error_message('Failed to mask '. $self->mask_reference_transcripts .' reference transcripts!');
        return;
    }
    return $tmp_bed_file;
};

1;
