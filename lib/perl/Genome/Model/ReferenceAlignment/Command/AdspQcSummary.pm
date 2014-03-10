package Genome::Model::ReferenceAlignment::Command::AdspQcSummary;

use strict;
use warnings;

use Genome;
use Spreadsheet::WriteExcel;

class Genome::Model::ReferenceAlignment::Command::AdspQcSummary {
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
        'Average Coverage',
        'Total Bases',
        'Total Mapped Bases',
        'Total Unique Mapped Bases',
        'Aligned %',
        'Unique %',
        'Error Rate',
        'Total # Reads',
        'Targeted Insert Length',
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
    my $flagstat_path = $result->bam_flagstat_path;
    my $flagstat = Genome::Model::Tools::Sam::Flagstat->parse_file_into_hashref($flagstat_path);
    my @instrument_data = $build->instrument_data;
    my @per_lane_alignments = $result->collect_individual_alignments;
    my ($total_bases, $total_mapped_bases, $total_unique_mapped_bases) = (0,0,0);
    my $mismatches = 0;
    my @inserts;
    my %readlengths;

    for my $lane (@per_lane_alignments) {
        $total_bases += $lane->total_base_count;
        $total_mapped_bases += $lane->total_aligned_base_count;
        $mismatches += $lane->total_base_count * $lane->instrument_data->filt_error_rate_avg / 100;
        push @inserts, $lane->instrument_data->library->original_insert_size;
        $readlengths{$lane->instrument_data->read_length} = 1;
    }

    if( scalar(keys %readlengths) > 1 ) {
        die "Multiple read lengths generated. NUMBERS WILL BE WRONG. I AM DYING TO LET YOU FIGURE IT OUT.\n";
    }

    #this is approximate. Running our c alignment stat tool to double check...
    $total_unique_mapped_bases = $total_mapped_bases - $flagstat->{reads_marked_duplicates} * (keys %readlengths)[0];

    return [
        $build->subject_name,
        scalar(@per_lane_alignments),
        $build_metrics{'haploid_coverage'},
        $total_bases,
        $total_mapped_bases,
        $total_unique_mapped_bases,
        sprintf("%0.02f",$flagstat->{reads_mapped}/$flagstat->{total_reads} * 100),
        sprintf("%0.02f",($flagstat->{reads_mapped} - $flagstat->{reads_marked_duplicates})/$flagstat->{total_reads} * 100),
        sprintf("%0.02f",$mismatches / $total_bases * 100),
        $flagstat->{total_reads},
        join(",", sort {$a <=> $b} @inserts)
        ];
}


1;

