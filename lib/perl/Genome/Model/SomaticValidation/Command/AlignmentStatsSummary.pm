package Genome::Model::SomaticValidation::Command::AlignmentStatsSummary;

use strict;
use warnings;

use Genome;

class Genome::Model::SomaticValidation::Command::AlignmentStatsSummary {
    is => 'Genome::SoftwareResult::StageableSimple',
    doc => 'Generate a spreadsheet, tsv file, of alignment metrics for Somatic Validation models.  Duplicate samples using the same processing profile will only be reported once.',
    has_input => [
        builds => {
            is => 'Genome::Model::Build::SomaticValidation',
            is_many => 1,
        },
    ],
    has_param => [
        haploid_coverage => {
            is => 'Boolean',
        },
        targeted_insert_length => {
            is => 'Boolean',
        },
    ],
    has_transient_optional => [
        _writer => {
            is => 'Genome::Utility::IO::SeparatedValueWriter',
        },
    ],
};

sub output_tsv_file_name {
    return "alignment_stats.tsv";
}

sub output_tsv_file_path {
    my $self = shift;
    return File::Spec->join($self->output_dir, $self->output_tsv_file_name);
}

sub _temp_file_path {
    my $self = shift;
    return File::Spec->join($self->temp_staging_directory, $self->output_tsv_file_name);
}

sub _run {
    my $self = shift;

    $self->_load_writer;

    # TODO: This should probably use the merged_alignment_result to make unique entries per sample and avoid duplicates.
    my %sample_to_pp;
    for my $build ($self->builds) {
        my $tumor_sample = $build->tumor_sample;
        if ($sample_to_pp{$tumor_sample->name}) {
            if ($sample_to_pp{$tumor_sample->name} eq $build->processing_profile->id) {
                $self->status_message('Duplicate tumor sample \''. $tumor_sample->name .'\' allowed, but only reported once in output.');
            } else {
                die('Duplicate tumor sample \''. $tumor_sample->name  .'\' found with different processing profile IDs.  Exiting!');
            }
        } else {
            $sample_to_pp{$tumor_sample->name} = $build->processing_profile->id;

            my $tumor_alignment_result = $build->merged_alignment_result;
            unless ($self->write_sample_metrics($tumor_sample,$tumor_alignment_result)) {
                die('Failed to write alignment stats for tumor sample: '. $tumor_sample->name);
            }
        }

        my $normal_sample = $build->normal_sample;
        if ($sample_to_pp{$normal_sample->name}) {
            if ($sample_to_pp{$normal_sample->name} eq $build->processing_profile->id) {
                $self->status_message('Duplicate normal sample \''. $normal_sample->name .'\' allowed, but only reported once in output.');
            } else {
                die('Duplicate normal sample \''. $normal_sample->name  .'\' found with different processing profile IDs.  Exiting!');
            }
        } else {
            $sample_to_pp{$normal_sample->name} = $build->processing_profile->id;
            my $normal_alignment_result = $build->control_merged_alignment_result;
            unless ($self->write_sample_metrics($normal_sample,$normal_alignment_result)) {
                die('Failed to write alignment stats for normal sample: '. $normal_sample->name);
            }
        }
    }
    $self->_writer->output->close;
    return 1;
}

sub _load_writer {
    my $self = shift;

    my @headers = ('subject_name',
                   '# of lanes (or "sequence events")',
                   'Total Bases',
                   'Total Mapped Bases',
                   'Total Unique Mapped Bases',
                   'Aligned %',
                   'Unique %',
                   'Error Rate',
                   'Total # Reads',
                   '%pairs mapping across chromosomes',
               );
    if ($self->haploid_coverage) {
        push @headers, 'Average Coverage';
    }
    if ($self->targeted_insert_length) {
        push @headers, 'Targeted Insert Length';
    }
    my $writer = Genome::Utility::IO::SeparatedValueWriter->create(
        output => $self->_temp_file_path,
        separator => "\t",
        headers => \@headers,
    );
    $self->_writer($writer);
    return 1;
}

sub write_sample_metrics {
    my $self = shift;
    my ($sample,$alignment_result) = @_;

    my $writer = $self->_writer;
    my $data = $self->_alignment_metrics_from_result($alignment_result);
    $data->{subject_name} = $sample->name;
    $writer->write_one($data);
    return 1;
}

sub _alignment_metrics_from_result {
    my $self = shift;
    my $result = shift;
    
    my $flagstat_path = $result->bam_flagstat_path;
    my $flagstat = Genome::Model::Tools::Sam::Flagstat->parse_file_into_hashref($flagstat_path);
    
    my @per_lane_alignments = $result->collect_individual_alignments;
    my ($total_bases, $total_mapped_bases, $total_unique_mapped_bases) = (0,0,0);
    my $mismatches = 0;
    my $haploid_coverage = 0;
    my @inserts;
    my %readlengths;

    for my $lane (@per_lane_alignments) {
        $total_bases += $lane->total_base_count;
        $total_mapped_bases += $lane->total_aligned_base_count;
        $mismatches += $lane->total_base_count * $lane->instrument_data->filt_error_rate_avg / 100;
        push @inserts, $lane->instrument_data->library->original_insert_size;
        $readlengths{$lane->instrument_data->read_length} = 1;
    }
    if ($self->haploid_coverage && ($total_mapped_bases > 0) ) {
        # Stolen from Genome::Model::ReferenceAlignment::Report::Mapcheck;
        my $coverage = Genome::Model::Tools::Sam::Coverage->create(
            aligned_reads_file=> $result->bam_path,
            reference_file => $result->reference_build->full_consensus_path('fa'),
            return_output => 1,
            coverage_command => $ENV{GENOME_SW} . '/samtools/bamcheck/bamcheck-v0.13/bam-check -q 1',
        );
        my $bam_coverage_report = $coverage->execute;
        if (defined($bam_coverage_report) ) {
            $self->debug_message("Bam coverage report successfully generated.");
            $self->debug_message("Bam coverage report string: \n".$bam_coverage_report);
        }  else {
            $self->error_message("Could not generate Bam coverage report.");
            die($self->error_message);
        }
        # Stolen from Genome::Model::ReferenceAlignment::Report::Summary;
        if ( ($bam_coverage_report =~ m/Average depth across all non-gap regions: (\S+)/g ) ||
                 ($bam_coverage_report =~ m/\nAverage Coverage:(\S+)/g ) ) {
            $haploid_coverage = $1 if defined($1);
        } else {
            die('Failed to parse haploid coverage from : '.  $bam_coverage_report);
        }
    }

    if( scalar(keys %readlengths) > 1 ) {
        die "Multiple read lengths generated. NUMBERS WILL BE WRONG. I AM DYING TO LET YOU FIGURE IT OUT.\n";
    }

    #this is approximate. Running our c alignment stat tool to double check...
    $total_unique_mapped_bases = $total_mapped_bases - $flagstat->{reads_marked_duplicates} * (keys %readlengths)[0];
    my $data = {
        '# of lanes (or "sequence events")' => scalar(@per_lane_alignments),
        'Total Bases' => $total_bases,
        'Total Mapped Bases' => $total_mapped_bases,
        'Total Unique Mapped Bases' => $total_unique_mapped_bases,
        'Aligned %' => sprintf("%0.02f",$flagstat->{reads_mapped}/$flagstat->{total_reads} * 100),
        'Unique %' => sprintf("%0.02f",($flagstat->{reads_mapped} - $flagstat->{reads_marked_duplicates})/$flagstat->{total_reads} * 100),
        'Error Rate' => sprintf("%0.02f",$mismatches / $total_bases * 100),
        'Total # Reads' => $flagstat->{total_reads},
        '%pairs mapping across chromosomes' => sprintf("%0.02f", ($flagstat->{'reads_mapped_in_interchromosomal_pairs'} / $flagstat->{'reads_mapped_in_pair'}) * 100 ),
    };
    if ($self->haploid_coverage) {
        $data->{'Average Coverage'} = $haploid_coverage;
    }
    if ($self->targeted_insert_length) {
        $data->{'Targeted Insert Length'} = join(",", sort {$a <=> $b} @inserts),
    }
    return $data;
}

1;
