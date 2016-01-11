package Genome::InstrumentData::AlignmentResult::Command::PicardRnaSeqMetrics;

use strict;
use warnings;

use Genome;
use Sys::Hostname;
use File::Path;
use Params::Validate ':types';

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
            doc => 'The version of Picard to use.',
        },
        picard_strand_specificity => {
            is => 'Text',
            doc => 'The transcript strand specificity used by Picard.',
            valid_values => Genome::Model::Tools::Picard::CollectRnaSeqMetrics->__meta__->property("strand_specificity")->valid_values,
            is_optional => 1,
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

sub required_files_from_annotation_build {
    (qw/ rRNA_MT_file annotation_file /);
}

sub verify_annotation_build_has_required_files {
    my ($self, $annotation_build, $reference_build) = Params::Validate::validate_pos(
        @_, {isa => __PACKAGE__}, {type => OBJECT}, {type => OBJECT},
    );

    my @missing_files;
    for my $file_method ( required_files_from_annotation_build() ) {
        push @missing_files, $file_method if not $self->file_from_annotation_build(
            $annotation_build, $reference_build, $file_method,
        );
    }

    if ( @missing_files ) {
        die $self->error_message("Cannot proceed! Missing required files from annotation build: @missing_files",);
    }

    return 1;
}

sub file_from_annotation_build {
    my ($self, $annotation_build, $reference_build, $file_method) = Params::Validate::validate_pos(
        @_, {isa => __PACKAGE__}, {type => OBJECT}, {type => OBJECT}, {type => SCALAR},
    );

    my $file = $annotation_build->$file_method('gtf', $reference_build->id ,0);
    return $file if $file and -s $file;

    $file = $annotation_build->$file_method('gtf', undef, 0);
    return $file if $file and -s $file;

    return;
}

sub execute {
    my $self = shift;

    $self->do_input_builds_have_required_files($self->annotation_build, $self->reference_build);

    my $metrics_directory = $self->metrics_directory;
    unless (-d $metrics_directory) {
        Genome::Sys->create_directory($metrics_directory);
    }

    my $annotation_build = $self->annotation_build;
    my $reference_build = $self->reference_build;
    my $alignment_result = $self->alignment_result;
    my $picard_version = $self->picard_version;

    my $bam_file = $alignment_result->get_bam_file;

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
    my $rRNA_MT_gtf_file = $self->file_from_annotation_build($annotation_build, $reference_build, 'rRNA_MT_file');
    my $seqdict_file = $reference_build->get_sequence_dictionary('sam',$reference_build->species_name,$picard_version);
    my $to_intervals_cmd = Genome::Model::Tools::Gtf::ToIntervals->execute(
        gtf_file => $rRNA_MT_gtf_file,
        seqdict_file => $seqdict_file,
        interval_file => $ribo_intervals,
    );
    unless ($to_intervals_cmd && $to_intervals_cmd->result) {
        $self->error_message('Failed to convert the rRNA_MT GTF file '. $rRNA_MT_gtf_file .' to intervals: '. $ribo_intervals);
        return;
    }

    # Get all mRNA annotation in RefFlat format
    my $mRNA_gtf_file = $self->file_from_annotation_build($annotation_build, $reference_build, 'annotation_file');
    my $to_ref_flat_cmd = Genome::Model::Tools::Gtf::ToRefFlat->execute(
        input_gtf_file => $mRNA_gtf_file,
        output_file => $mRNA_ref_flat,
    );
    unless ($to_ref_flat_cmd && $to_ref_flat_cmd->result) {
        $self->error_message('Failed to convert the all_sequences GTF file '. $mRNA_gtf_file .'  to RefFlat: '. $mRNA_ref_flat);
        return;
    }

    # This is wasteful, but required since BAM sort order does not match our FASTA chromosome order
    my $tmp_bam_file = Genome::Sys->create_temp_file_path();
    my $reorder_sam_cmd = Genome::Model::Tools::Picard::ReorderSam->execute(
        input_file => $bam_file,
        output_file => $tmp_bam_file,
        reference_file => $reference_fasta_file,
        use_version => $picard_version,
    );
    unless ($reorder_sam_cmd && $reorder_sam_cmd->result) {
        $self->error_message('Failed to reorder BAM file: '. $bam_file);
        return;
    }
    my %picard_params = (
        input_file => $tmp_bam_file,
        output_file => $metrics_output_file,
        refseq_file => $reference_fasta_file,
        ribosomal_intervals_file => $ribo_intervals,
        ref_flat_file => $mRNA_ref_flat,
        use_version => $picard_version,
        chart_output => $chart_output_file,
    );
    if ($self->picard_strand_specificity) {
        $picard_params{strand_specificity} = $self->picard_strand_specificity;
    }
    my $collect_rna_seq_metrics_cmd = Genome::Model::Tools::Picard::CollectRnaSeqMetrics->execute(%picard_params);
    unless ($collect_rna_seq_metrics_cmd && $collect_rna_seq_metrics_cmd->result) {
        $self->error_message('Failed to run Picard CollectRnaSeqMetrics for alignment result: '. $alignment_result->id);
        return;
    }
    my $plot_rna_seq_metrics_cmd = Genome::Model::Tools::Picard::PlotRnaSeqMetrics->execute(
        input_file => $metrics_output_file,
        output_file => $pie_chart_file,
        label => $alignment_result->id,
    );
    unless ($plot_rna_seq_metrics_cmd && $plot_rna_seq_metrics_cmd->result) {
        $self->error_message('Failed to run PlotRnaSeqMetrics for alignment result: '. $alignment_result->id);
        return;
    }
    return 1;
}


1;
