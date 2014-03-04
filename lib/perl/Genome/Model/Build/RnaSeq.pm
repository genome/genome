package Genome::Model::Build::RnaSeq;

use strict;
use warnings;

use Genome;
use File::Path 'rmtree';

my $DEFAULT_RIBO_INTERVALS_FILENAME = 'Picard_ribo.intervals';
my $DEFAULT_MRNA_REF_FLAT_FILENAME = 'Picard_mRNA.refFlat';
my $DEFAULT_METRICS_OUTPUT_FILENAME = 'PicardRnaSeqMetrics.txt';
my $DEFAULT_CHART_OUTPUT_FILENAME = 'PicardRnaSeqChart.pdf';
my $DEFAULT_PIE_CHART_FILENAME = 'PicardRnaSeqMetrics.png';

class Genome::Model::Build::RnaSeq {
    is => 'Genome::Model::Build',
    has => [
        annotation_build => {
            is => "Genome::Model::Build::ImportedAnnotation",
            is_input => 1,
            is_optional => 1,
        },
        reference_sequence_build => {
            is => "Genome::Model::Build::ReferenceSequence",
            is_input => 1
        },
    ]
};

sub accumulated_alignments_directory {
    my $self = shift;
    return $self->data_directory . '/alignments';
}

sub alignment_stats_directory {
    my $self = shift;
    # this is really the name of the symlink to the software result
    return $self->data_directory . '/alignment_stats';
}

sub alignment_stats_file {
    my $self = shift;
    # TopHat1
    if ($self->processing_profile->read_aligner_version =~ /^1/) {
        return $self->accumulated_alignments_directory .'/alignment_stats.txt';
    } else {
        if (-e $self->alignment_stats_directory) {
            # A new directory was created since picard metrics will be moved to a software result and symlinked to the build
            return $self->alignment_stats_directory .'/alignment_stats.txt';
        } else {
            # For a brief one-week period the file was linked to the metrics directory
            # This was before picard metrics were a software result
            return  $self->metrics_directory .'/alignment_stats.txt';
        }
    }

}

sub coverage_directory {
    my $self = shift;
    return $self->data_directory . '/coverage';
}

sub metrics_directory {
    my $self = shift;
    return $self->data_directory . '/metrics';
}

sub junctions_directory {
    my $self = shift;
    return $self->data_directory . '/junctions';
}

sub default_ribo_intervals_filename {
    my $class = shift;
    return $DEFAULT_RIBO_INTERVALS_FILENAME;
}

sub picard_rna_seq_ribo_intervals {
    my $self = shift;
    return $self->metrics_directory .'/'. $self->picard_rna_seq_ribo_intervals;
}

sub default_mRNA_ref_flat_filename {
    my $class = shift;
    return $DEFAULT_MRNA_REF_FLAT_FILENAME;
}

sub picard_rna_seq_mRNA_ref_flat {
    my $self = shift;
    return $self->metrics_directory .'/'. $self->default_mRNA_ref_flat_filename;
}

sub default_metrics_output_filename {
    my $class = shift;
    return $DEFAULT_METRICS_OUTPUT_FILENAME;
}

sub picard_rna_seq_metrics {
    my $self = shift;
    return $self->metrics_directory .'/'. $self->default_metrics_output_filename;
}

sub default_chart_output_filename {
    my $class = shift;
    return $DEFAULT_CHART_OUTPUT_FILENAME;
}

sub picard_rna_seq_chart {
    my $self = shift;
    return $self->metrics_directory .'/'. $self->default_chart_output_filename;
}

sub default_pie_chart_filename {
    my $class = shift;
    return $DEFAULT_PIE_CHART_FILENAME;
}

sub picard_rna_seq_pie_chart {
    my $self = shift;
    return $self->metrics_directory .'/'. $self->picard_rna_seq_pie_chart;
}

sub accumulated_alignments_disk_allocation {
    my $self = shift;

    my $align_event = Genome::Model::Event::Build::RnaSeq::AlignReads->get(
        model_id=>$self->model->id,
        build_id=>$self->build_id
    );

    return if (!$align_event);

    my $disk_allocation = Genome::Disk::Allocation->get(owner_class_name=>ref($align_event), owner_id=>$align_event->id);

    return $disk_allocation;
}

sub accumulated_fastq_directory {
    my $self = shift;
    return $self->data_directory . '/fastq';
}

sub expression_directory {
    my $self = shift;
    my $mode = shift;

    $mode = 'expression' unless $mode;

    if ($mode eq 'expression') {
        return $self->accumulated_expression_directory;
    } elsif ($mode eq 'de novo' or $mode eq 'de_novo') {
        return $self->de_novo_directory;
    } elsif ($mode eq 'reference guided' or $mode eq 'reference_guided') {
        return $self->reference_guided_directory;
    } elsif ($mode eq 'reference only' or $mode eq 'reference_only') {
        return $self->reference_only_directory;
    } else {
        die $self->error_message("'$mode' is not a valid mode for annotation_reference_transcripts_mode.");
    }
}

sub accumulated_expression_directory {
    my $self = shift;
    return $self->data_directory . '/expression';
}

sub de_novo_directory {
    my $self = shift;
    return $self->data_directory . '/de_novo';
}

sub reference_guided_directory {
    my $self = shift;
    return $self->data_directory . '/reference_guided';
}

sub reference_only_directory {
    my $self = shift;
    return $self->data_directory . '/reference_only';
}

sub merged_alignment_result {
    my $self = shift;
    return $self->alignment_result;
}

sub alignment_result {
    my $self = shift;

    my @u = Genome::SoftwareResult::User->get(user_id => $self->build_id);
    my $alignment_class = "";
    if($self->processing_profile->read_aligner_name eq 'tophat' && $self->processing_profile->read_aligner_version =~ /^1/){
      $alignment_class = Genome::InstrumentData::AlignmentResult->_resolve_subclass_name_for_aligner_name($self->processing_profile->read_aligner_name);
    }else{
      $alignment_class = "Merged";
    }
    my $alignment = join('::', 'Genome::InstrumentData::AlignmentResult', $alignment_class)->get([map($_->software_result_id, @u)]);
    return $alignment;
}

sub alignment_result_with_lock {
    my $self = shift;

    return $self->_fetch_alignment_result('get_with_lock');
}

sub generate_alignment_result {
    my $self = shift;

    return $self->_fetch_alignment_result('get_or_create');
}

sub _fetch_alignment_result {
    my $self = shift;
    my $mode = shift;

    my @instrument_data_inputs = $self->instrument_data_inputs;
    my ($params) = $self->model->params_for_alignment(@instrument_data_inputs);

    my $alignment_class = Genome::InstrumentData::AlignmentResult->_resolve_subclass_name_for_aligner_name($self->model->read_aligner_name);
    my $alignment = join('::', 'Genome::InstrumentData::AlignmentResult', $alignment_class)->$mode(
        %$params,
    );

    return $alignment;
}

sub delete {
    my $self = shift;
    
    # if we have an alignments directory, nuke it first since it has its own allocation
    if (-e $self->accumulated_alignments_directory ||
        -e $self->accumulated_fastq_directory ||
        -e $self->accumulated_expression_directory ||
        -e $self->de_novo_directory ||
        -e $self->reference_guided_directory ||
        -e $self->reference_only_directory
    ) {
        unless($self->eviscerate()) {
            my $eviscerate_error = $self->error_message();
            $self->error_message("Eviscerate failed: $eviscerate_error");
            return;
        };
    }
    
    $self->SUPER::delete(@_);
}

# nuke the accumulated alignment directory
sub eviscerate {
    my $self = shift;
    
    $self->debug_message('Entering eviscerate for build:' . $self->id);


    if($self->alignment_result) {
        my $alignment_result = $self->alignment_result;

        if (-l $self->accumulated_alignments_directory && readlink($self->accumulated_alignments_directory) eq $alignment_result->output_dir) {
           $self->debug_message("Unlinking symlink to alignment result: " . $self->accumulated_alignments_directory);
            unless(unlink($self->accumulated_alignments_directory)) {
                $self->error_message("could not remove symlink to alignment result path");
                return;
            }
        }

        my @users = $alignment_result->users(user => $self);
        map($_->delete, @users);
        $self->debug_message('Removed self as user of alignment result.');
    } else {
        my $alignment_alloc = $self->accumulated_alignments_disk_allocation;
        my $alignment_path = ($alignment_alloc ? $alignment_alloc->absolute_path :  $self->accumulated_alignments_directory);

        if (!-d $alignment_path && !-l $self->accumulated_alignments_directory) {
            $self->debug_message("Nothing to do, alignment path doesn't exist and this build has no alignments symlink.");
        }

        $self->debug_message("Removing tree $alignment_path");
        if (-d $alignment_path) {
            rmtree($alignment_path);
            if (-d $alignment_path) {
                $self->error_message("alignment path $alignment_path still exists after evisceration attempt, something went wrong.");
                return;
            }
        }

        if ($alignment_alloc) {
            unless ($alignment_alloc->deallocate) {
                $self->error_message("could not deallocate the alignment allocation.");
                return;
            }
        }

        if (-l $self->accumulated_alignments_directory && readlink($self->accumulated_alignments_directory) eq $alignment_path ) {
            $self->debug_message("Unlinking symlink: " . $self->accumulated_alignments_directory);
            unless(unlink($self->accumulated_alignments_directory)) {
                $self->error_message("could not remove symlink to deallocated accumulated alignments path");
                return;
            }
        }
    }

    my $fastq_directory            = $self->accumulated_fastq_directory;
    my $expression_directory       = $self->accumulated_expression_directory;
    my $de_novo_directory          = $self->de_novo_directory;
    my $reference_guided_directory = $self->reference_guided_directory;
    my $reference_only_directory   = $self->reference_only_directory;

    my %directories = (
        fastq            => $fastq_directory,
        expression       => $expression_directory,
        de_novo          => $de_novo_directory,
        reference_guided => $reference_guided_directory,
        reference_only   => $reference_only_directory,
    );

    for my $item (keys %directories) {
        my $directory = $directories{$item};
        if (-d $directory) {
            $self->debug_message("removing $item directory");
            rmtree($directory);
            if (-d $directory) {
                $self->error_message("$item path $directory still exists after evisceration attempt, something went wrong.");
                return;
            }
        }
    }

    return 1;
}

sub workflow_name {
    my $self = shift;
    return $self->build_id;
}

sub workflow_instances {
    my $self = shift;

    my @instances = $self->SUPER::workflow_instances;

    unless(@instances) {
        @instances = Workflow::Operation::Instance->get(
            name => $self->SUPER::workflow_name, #older profiles were staged
        );
    }
    return @instances;
}

sub ensure_annotation_build_provided {
    my $self = shift;
    my @tags = ();
    unless ( ($self->model->annotation_reference_transcripts_mode eq 'de novo') or $self->annotation_build ) {
        push @tags, UR::Object::Tag->create(
            type => 'invalid',
            properties => [qw/ annotation_build /],
            desc => "Processing Profile calls for " . $self->model->annotation_reference_transcripts_mode . " mode, but this model does not have an annotation_build set",
        );
    }
    return @tags;
}

sub validate_for_start_methods {
    my $self = shift;
    my @methods = $self->SUPER::validate_for_start_methods;
    push @methods, 'ensure_annotation_build_provided';
    return @methods;
}

sub reference_being_replaced_for_input {
    my $self = shift;
    my $input = shift;

    return unless $input;
    return if $self->processing_profile->read_aligner_name eq 'imported';

    #we're going to realign, so any existing reference on the data is unimportant
    return 1 if $input->name eq 'instrument_data';

    return;
}

# DIFFING RELATED FUNCTIONS

sub regex_for_custom_diff {
    my $self = shift;

    my @regexes_from_base_class = $self->SUPER::regex_for_custom_diff;
    return (@regexes_from_base_class,
        'bam_without_regard_to_header' => '\.bam$',
        'via_md5' => '\.fa$',
        'via_md5' => '\.ebwt$',
    );
}

sub diff_via_md5 {
    my ($self, $first_file, $second_file) = @_;

    my $first_md5  = `md5sum $first_file`;
    my $second_md5 = `md5sum $second_file`;
    return 1 if $first_md5 eq $second_md5;
    return 0;
}

sub diff_bam_without_regard_to_header {
    my ($self, $first_file, $second_file) = @_;

    my $first_md5  = `samtools view $first_file | md5sum`;
    my $second_md5 = `samtools view $second_file | md5sum`;
    return 1 if $first_md5 eq $second_md5;
    return 0;
}

sub files_ignored_by_diff {
    my $self = shift;

    return (
        '.fai$',
        'build.xml',
        '.pdf$',
        '.png$',
        '.log$',
        '.metrics$',
        '.html$',
        '.zip$',
        '.bai',
        '.swp',
        'alignments/[[:xdigit:]]+(?:_merged_rmdup)?.bam.md5$',
        'R.stderr$',
        '_fastqc/fastqc_data.txt',
        '_fastqc/summary.txt',
        '-PicardGC_metrics.txt',
        '-PicardGC_summary.txt',
        File::Spec->join('metrics', 'PicardRnaSeqMetrics.txt'),

        'fusions.*Index.*',
        'fusions.*runconfig.xml',
        'fusions.*sorted_aligned_reads.bam.*',
    );
}

sub dirs_ignored_by_diff {
    my $self = shift;

    my @directories = qw(
        logs
        reports
        reference_guided
        reference_only
        de_novo
        expression
    );
    return @directories;
}

sub regex_files_for_diff {
    my $self = shift;

    return ($self->alignment_regexes, $self->junctions_regexes, $self->bam_qc_regexes);

}
sub alignment_regexes {
    my $self = shift;
    return qw(
        alignments/[[:xdigit:]]+(?:_merged_rmdup)?.bam$
        alignments/[[:xdigit:]]+(?:_merged_rmdup)?.bam.flagstat$
    );
}

sub junctions_regexes {
    my $self = shift;

    return qw(
        junctions/AlignmentResult_[[:xdigit:]]+_junctions.bed$
    );
}

sub bam_qc_regexes {
    my $self = shift;
    my @files = qw(
        -AlignmentLengthDistribution.tsv
        -ErrorRate.tsv
        -ReadLengthDistribution.tsv
        -ReadLengthSummary.tsv
        .bam
        .bam.bai
        .bam.flagstat
    );

    my @regexes = map {sprintf('bam-qc/[[:xdigit:]]+%s$', $_)} @files;
    return @regexes;
}


1;

