package Genome::Model::Tools::Capture::EvaluateExomeBreadth;

use strict;
use warnings;

my $DEFAULT_DELIMITER = ':';
my $DEFAULT_ANNOTATION_FILE = '/gscmnt/sata132/techd/solexa/jwalker/RNAseq/annotation/hg18/NCBI-human.combined-annotation-54_36p_v2/NCBI-human.combined-annotation_54_36p_v2.bed';

class Genome::Model::Tools::Capture::EvaluateExomeBreadth {
    is => ['Command'],
    has => [
        include_list => {
            is => 'Text',
            doc => 'A list of gene names to include.',
            is_optional => 1,
        },
        exome_bed_file => {
            is => 'Text',
            doc => 'The Exome product BED format file.  No "chr"',
        },
        output_file => {
            is => 'Text',
            doc => 'The output file.',
        },
        annotation_bed_file => {
            is => 'Text',
            is_optional => 1,
            doc => 'The annotation BED format file.  With features named like $GENE:$TRANSCRIPT:$TYPE:$ORDINAL.',
            default_value => $DEFAULT_ANNOTATION_FILE,
        },
        coverage_of => {
            is => 'Text',
            is_optional => 1,
            doc => 'A comma delimited list of transcript sub structures to evaluate target breadth.  Valid values: utr_exon,cds_exon,intron,rna',
            default_value => 'cds_exon',
        },
        report_by => {
            is => 'Text',
            is_optional => 1,
            default_value => 'gene',
            valid_values => ['exome','gene','transcript','exon'],
        },
        delimiter => {
            is => 'Text',
            doc => 'The character that delimits GENE, TRANSCRIPT, TYPE, and ORDINAL',
            is_optional => 1,
            default_value => $DEFAULT_DELIMITER,
            valid_values => [':','_','.'],
        },
        bedtools_version => {
            is => 'Text',
            is_optional => 1,
            doc => 'The version of BEDTools to use.',
            default_value => Genome::Model::Tools::BedTools->default_bedtools_version,
        },
        samtools_version => {
            is => 'Text',
            is_optional => 1,
            doc => 'The version of samtools to use.',
            default_value => Genome::Model::Tools::Sam->default_samtools_version,
        },
        minimum_breadth_filters => {
            is => 'Text',
            is_optional => 1,
            doc => 'A comma-delimited list of minimum breadths as percentages.',
            default_value => '100,90,80',
        },
    ],
};

sub default_delimiter {
    return $DEFAULT_DELIMITER;
}

sub execute {
    my $self = shift;

    my ($exome_basename,$exome_dirname,$exome_suffix) = File::Basename::fileparse($self->exome_bed_file,qw/\.bed/);
    unless ($exome_basename && $exome_suffix) {
        die('Failed to parse BED file path:  '. $self->exome_bed_file);
    }
    my ($annotation_basename,$annotation_dirname,$annotation_suffix) = File::Basename::fileparse($self->annotation_bed_file,qw/\.bed/);
    unless ($annotation_basename && $annotation_suffix) {
        die('Failed to parse BED file path:  '. $self->annotation_bed_file);
    }
    my $sorted_bam_file = Genome::Sys->create_temp_file_path($exome_basename .'_sorted.bam');
    my $BedToSortedBam = Genome::Model::Tools::BedTools::BedToSortedBam->execute(
        input_file => $self->exome_bed_file,
        output_file => $sorted_bam_file,
        use_version => $self->bedtools_version,
        samtools_version => $self->samtools_version,
    );
    unless ($BedToSortedBam) {
        die('Failed to convert BED file to sorted BAM file!');
    }
    my $limited_bed_file = Genome::Sys->create_temp_file_path($annotation_basename .'_limited_to_list.bed');
    unless (Genome::Model::Tools::Bed::Limit->execute(
        gene_list => $self->include_list,
        input_bed_file => $self->annotation_bed_file,
        output_bed_file => $limited_bed_file,
        feature_types => $self->coverage_of,
    )) {
        die('Failed to limit BED file '. $self->annotation_bed_file .' to the genes in '. $self->include_list);
    }
    my $merged_bed_file = Genome::Sys->create_temp_file_path($annotation_basename .'_merged_by_'. $self->report_by .'.bed');
    unless ($self->report_by eq 'exon') {
        my $MergeBy = Genome::Model::Tools::BedTools::MergeBy->create(
            merge_by => $self->report_by,
            input_file => $limited_bed_file,
            output_file => $merged_bed_file,
            delimiter => $self->delimiter,
            use_version => $self->bedtools_version,
        );
        unless ($MergeBy->execute) {
            die('Failed to merge annotation by '. $self->report_by);
        }
    } else {
        $merged_bed_file = $limited_bed_file;
    }
    my $stats_file = Genome::Sys->create_temp_file_path($annotation_basename .'_merged_by_'. $self->report_by .'_STATS.tsv');;
    my $RefCov = Genome::Model::Tools::BioSamtools::RefCov->execute(
        bed_file => $merged_bed_file,
        bam_file => $sorted_bam_file,
        stats_file => $stats_file,
    );
    unless ($RefCov) {
        die('Failed to genarate coverage!');
    }
    my $MergeStats = Genome::Model::Tools::BioSamtools::MergeStats->execute(
        stats_file => $stats_file,
        output_file => $self->output_file,
        merge_by => $self->report_by,
        delimiter => $self->delimiter,
        minimum_breadth_filters => $self->minimum_breadth_filters,
    );
    unless ($MergeStats) {
        die ('Failed to merge stats by '. $self->report_by);
    }
    return 1;
}


1;
