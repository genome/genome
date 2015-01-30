package Genome::Model::ReferenceAlignment::Command::AdspExomeQcSummary;

use strict;
use warnings;

use Genome;
use Spreadsheet::WriteExcel;
use Data::Dumper;

class Genome::Model::ReferenceAlignment::Command::AdspExomeQcSummary {
    is => 'Command::V2',
    doc => "Generate ADSP requested QC metrics as an Excel Spreadsheet.",
    has => [
        builds => {
            is => 'Genome::Model::Build::ReferenceAlignment',
            shell_args_position => 1,
            is_many => 1,
            doc => 'List of builds to report on',
        },
        output_file => {
            is => 'Text',
            is_optional => 1,
            doc => 'The file to write output into.'
        }
    ],
};


sub help_detail {
    return "Generate tab-delimited (*.tsv) files with BamQc metrics summarized by instrument data id for the provided build.";
}


sub execute {
    my $self  = shift;

    my $filename = $self->output_file;
    unless($filename =~ /.*\.xls$/) {
        $self->status_message("No xls extension detected. Appending to the filename.");
        $filename .= ".xls";
    }
    my $workbook = Spreadsheet::WriteExcel->new($filename);
    my $sheet = $workbook->add_worksheet();
    my @header = (
        'sample',
        '# of lanes (or "sequence events")',
        'Total Bases',
        'Total Mapped Bases',
        'Total Unique Mapped Bases',
        'Aligned %',
        'Unique %',
        'Error Rate',
        'Total # Reads',
        'Targeted Insert Length',
        'Percent Target Bases > 20X',
        'Total Unique On-Target Bases',
        'Average On-Target Depth',
    );
    my $row = 0;
    $sheet->write_row($row++,0, \@header);

    for my $build ($self->builds) {
        if($build->status ne 'Succeeded') {
            my $msg = sprintf("Filtering out build %s (model %s) because its status is %s, not succeeded.", $build->id, $build->model->name, $build->status);
            $self->status_message($msg);
        }
        else {
            $sheet->write_row($row++,0, $self->build_stats($build));
        }
    }
    

    return 1;
}

sub build_stats {
    my ($self,$build) = @_;
    my @metrics;
    my %build_metrics = map { $_->name => $_->value} $build->metrics;
    my $result = $build->merged_alignment_result;
    unless(defined $result) {
        warn "Unable to find a merged alignment result for " . $build->id ." (" . $build->model->name ."). Skipping this model and build. You should investigate further.\n";
        next;
    }
    my $flagstat_path = $result->bam_flagstat_path;
    my $flagstat = Genome::Model::Tools::Sam::Flagstat->parse_file_into_hashref($flagstat_path);
    my @instrument_data = $build->instrument_data;
    my @per_lane_alignments = $build->alignment_results;
    my ($total_bases, $total_mapped_bases, $total_unique_mapped_bases) = (0,0,0);
    my $mismatches = 0;
    my @inserts;

    for my $lane (@per_lane_alignments) {
        my $avg_mismatch_rate = ($lane->metric(metric_name => 'read_1_pct_mismatch')->metric_value + $lane->metric(metric_name => 'read_2_pct_mismatch')->metric_value)/200;
        $mismatches += $lane->total_base_count * $avg_mismatch_rate;
        push @inserts, $lane->metric(metric_name => 'median_insert_size')->metric_value;
    }

    my @target_stats;
    my $coverage_stats_summary_hash_ref = $build->coverage_stats_summary_hash_ref();
    my $wingspan_coverage_stats_summary_hash_ref = $coverage_stats_summary_hash_ref->{0};
    my $alignment_summary_hash_ref = $build->reference_coverage_result->alignment_summary_v2_hash_ref();
    # ADSP wants 20X
    my $data = $wingspan_coverage_stats_summary_hash_ref->{20};
    my $pct_covered_greater_than_20X = $wingspan_coverage_stats_summary_hash_ref->{20}->{pc_target_space_covered};
    my $total_unique_on_target = $alignment_summary_hash_ref->{0}->{unique_target_aligned_bp};
    $total_bases = $alignment_summary_hash_ref->{0}->{total_bp};
    $total_mapped_bases = $alignment_summary_hash_ref->{0}->{total_aligned_bp};
    $total_unique_mapped_bases = $total_mapped_bases - $alignment_summary_hash_ref->{0}->{total_duplicate_bp}; 
    my $average_on_target_depth = $total_unique_on_target / $wingspan_coverage_stats_summary_hash_ref->{20}->{target_base_pair};

    @target_stats = ($pct_covered_greater_than_20X, $total_unique_on_target, sprintf("%0.02f", $average_on_target_depth));

    my @stats = (
        $build->subject_name,
        scalar(@per_lane_alignments),
        $total_bases,
        $total_mapped_bases,
        $total_unique_mapped_bases,
        sprintf("%0.02f",$flagstat->{reads_mapped}/$flagstat->{total_reads} * 100),
        sprintf("%0.02f",($flagstat->{reads_mapped} - $flagstat->{reads_marked_duplicates})/$flagstat->{total_reads} * 100),
        sprintf("%0.02f",$mismatches / $total_bases * 100),
        $flagstat->{total_reads},
        join(",", sort {$a <=> $b} @inserts),
        @target_stats,
    );
    return \@stats;
}


1;

