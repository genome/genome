package Genome::InstrumentData::AlignmentResult::Merged::SpliceJunctionSummary;

use strict;
use warnings;

use Genome;
use Sys::Hostname;
use File::Path;

class Genome::InstrumentData::AlignmentResult::Merged::SpliceJunctionSummary {
    is => ['Genome::SoftwareResult::Stageable'],
    has_input => [
        alignment_result_id => {
            is => 'Number',
            doc => 'ID of the result for the alignment data upon which to run coverage stats',
        },
    ],
    has_param => [
        bedtools_version => {
            is => 'Text',
            doc => 'The version of Picard to use.',
        },
    ],
    has_metric => [
        _log_directory => {
            is => 'Text',
            doc => 'Path where workflow logs were written',
        },
        #many other metrics exist--see sub _generate_metrics
    ],
    has => [
        alignment_result => {
            is => 'Genome::InstrumentData::AlignmentResult::Merged',
            id_by => 'alignment_result_id',
            doc => 'the alignment data upon which to run coverage stats',
        },
    ],
    has_transient_optional => [
        log_directory => {
            is => 'Text',
            doc => 'Path to write logs from running the workflow',
        },
    ],
};

sub resolve_allocation_subdirectory {
    my $self = shift;

    my $hostname = hostname;
    my $user = $ENV{'USER'};
    my $base_dir = sprintf("splice_junction_summary-%s-%s-%s-%s", $hostname, $user, $$, $self->id);

    # TODO: the first subdir is actually specified by the disk management system.
    my $directory = join('/', 'build_merged_alignments','splice_junction_summary',$base_dir);
    return $directory;
}

sub resolve_allocation_disk_group_name {
    $ENV{GENOME_DISK_GROUP_MODELS};
}

sub _staging_disk_usage {
    #need the allocation created in advance for this process
    return 5_000_000; #TODO better estimate
}

sub _working_dir_prefix {
    return 'splice-junction-summary';
}

sub create {
    my $class = shift;
    my $self = $class->SUPER::create(@_);
    return unless ($self);

    $self->_prepare_staging_directory;

    # Reference inputs
    my $reference_build = $self->alignment_result->reference_build;
    my $reference_fasta_file = $reference_build->full_consensus_path('fa');
    die $self->error_message("Reference FASTA File ($reference_fasta_file) is missing") unless -s $reference_fasta_file;
    
    # Annotation inputs
    my $annotation_build = $self->alignment_result->annotation_build;
    my $annotation_gtf_file = $annotation_build->annotation_file('gtf',$reference_build->id);
    unless(-s $annotation_gtf_file) {
        $annotation_gtf_file = $annotation_build->annotation_file('gtf');
    }
    my $annotation_file_basename = $annotation_build->name;
    $annotation_file_basename =~ s/\//-/g;

    # Alignment inputs
    my @alignments = $self->alignment_result->collect_individual_alignments;
    $self->debug_message('Merging '. scalar(@alignments) .' per lane junctions BED12 files...');

    my @alignment_junctions_bed12_files;
    for my $alignment (@alignments) {
        my $alignment_junctions_bed12_file = $alignment->output_dir .'/junctions.bed';
        my $id = $alignment->id;
        Genome::Sys->create_symlink($alignment_junctions_bed12_file,$self->temp_staging_directory .'/AlignmentResult_'. $id .'_junctions.bed');
        push @alignment_junctions_bed12_files, $alignment_junctions_bed12_file;
    }

    my $merged_junctions_bed12_file = $self->temp_staging_directory  .'/junctions.bed';
    if (scalar(@alignment_junctions_bed12_files) == 1) {
        Genome::Sys->create_symlink($alignment_junctions_bed12_files[0], $merged_junctions_bed12_file);
    } else {
        unless (Genome::Model::Tools::Bed::MergeBed12Junctions->execute(
            input_files => \@alignment_junctions_bed12_files,
            output_file => $merged_junctions_bed12_file,
            bedtools_version => $self->bedtools_version,
        )) {
            die();
        }
    }

    if(Genome::Sys->line_count($merged_junctions_bed12_file) < 2) {
        $self->warning_message("No splice junctions found for alignments.  Skipping run.");
        $self->_finalize_allocation();
        return $self;
    }

    my $log_dir = $self->log_directory;
    unless($log_dir) {
        $log_dir = '' . $self->temp_staging_directory;
    }
    $self->_log_directory($log_dir);

    my %params = (
        output_directory => $self->temp_staging_directory,
        annotation_name => $annotation_file_basename,
        annotation_gtf_file => $annotation_gtf_file,
        reference_fasta_file => $reference_fasta_file,
        observed_junctions_bed12_file => $merged_junctions_bed12_file,
        bedtools_version => $self->bedtools_version,
    );

    my $cmd = Genome::Model::Tools::Transcriptome::SpliceJunctionSummary->create(%params);
    unless($cmd->execute) {
        die('Failed to run SpliceJunctionSummary tool with params: '. Data::Dumper::Dumper(%params));
    }

    $self->_finalize_allocation();

    $self->_generate_metrics();

    return $self;
}

sub _finalize_allocation {
    my $self = shift;

    $self->_prepare_output_directory;
    $self->_promote_data;
    $self->_reallocate_disk_allocation;

    return 1;
}

sub _generate_metrics {
    my $self = shift;
    $self->debug_message('Currently no metrics are saved for: '. __PACKAGE__);
    #my $metrics = shift;
    
    #for my $type_label (keys %{$metrics}) {
        # Currently, do not store insert size, quality by cycle, or mean quality histograms
    #    if ($type_label =~ /Histogram/) { next; }
        # Currently, do not store the GcBiasMetrics, 100 Windows of normalized coverage
     #   if ($type_label eq 'GcBiasMetrics') { next; }
        
      #  my $type_metrics = $metrics->{$type_label};
        # The FlagstatMetrics hashref are one-level
       # if ($type_label eq 'FlagstatMetrics') {
        #    for my $metric_label (keys %{$type_metrics}) {
         #       my $metric_key = sprintf('bam_qc-%s-%s',$type_label,$metric_label);
          #      $self->add_metric(metric_name => $metric_key, metric_value => $type_metrics->{$metric_label});
           # }
        #} else {
            # All other hashrefs are considered to habe two-levels, the first level of the hashref being the key on which lines of metrics are differentiated
         #   for my $key (keys %{$type_metrics}) {
          #      for my $metric_label (keys %{$type_metrics->{$key}}) {
          #          my $metric_key = sprintf('bam_qc-%s-%s-%s',$type_label,$key,$metric_label);
          #          $self->add_metric(metric_name => $metric_key, metric_value => $type_metrics->{$key}->{$metric_label});
           #     }
           # }
        #}
   # }

    return 1;
}


1;
