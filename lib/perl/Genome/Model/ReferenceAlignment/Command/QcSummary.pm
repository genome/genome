package Genome::Model::ReferenceAlignment::Command::QcSummary;

use strict;
use warnings;

use Genome;
use Spreadsheet::WriteExcel;

class Genome::Model::ReferenceAlignment::Command::QcSummary {
    is => 'Command::V2',
    doc => "Generate QC metrics as an Excel Spreadsheet.",
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
    return "Generate an Excel (*.xls) file with Qc metrics summarized for the provided builds.";
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
        'sample_name',
        'build_id',
        'instrument_data',
        'instrument_data_count',
        'bam_file',
        'total_reads',
        'total_bases',
        'total_aligned_reads',
        'total_duplicate_reads',
        'alignment_rate',
        'duplication_rate',
        'error_rate',
        'haploid_coverage',
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

    my %build_metrics = map { $_->name => $_->value} $build->metrics;
    
    my $result = $build->merged_alignment_result;

    my $flagstat_path = $result->bam_flagstat_path;
    my $flagstat = Genome::Model::Tools::Sam::Flagstat->parse_file_into_hashref($flagstat_path);
    my $total_duplicate_reads = $flagstat->{reads_marked_duplicates};
    
    my ($total_reads, $total_bases, $total_aligned_reads,$mismatches) = (0,0,0,0);
    my @instrument_data = $build->instrument_data;
    for my $instrument_data (@instrument_data) {
        my @instrument_data_results = $build->alignment_results_for_instrument_data($instrument_data);
        unless (@instrument_data_results == 1) { die;}
        my $instrument_data_result = $instrument_data_results[0];

        $total_reads += $instrument_data_result->total_read_count;
        $total_bases += $instrument_data_result->total_base_count;
        $total_aligned_reads += $instrument_data_result->total_aligned_read_count;
        $mismatches += $instrument_data_result->total_base_count * $instrument_data_result->instrument_data->filt_error_rate_avg / 100;
    }

    return [
        $build->model->subject->name,
        $build->id,
        join(',',map {$_->id} @instrument_data),
        scalar(@instrument_data),
        $result->bam_path,
        $total_reads,
        $total_bases,
        $total_aligned_reads,
        $total_duplicate_reads,
        sprintf("%0.04f",($total_aligned_reads / $total_reads)),
        sprintf("%0.04f",($total_duplicate_reads / $total_reads)),
        sprintf("%0.04f",($mismatches / $total_bases)),
        $build_metrics{'haploid_coverage'},
    ];
}


1;

