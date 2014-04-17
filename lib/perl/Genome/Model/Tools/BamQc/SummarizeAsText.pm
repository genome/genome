package Genome::Model::Tools::BamQc::SummarizeAsText;

use strict;
use warnings;

use Genome;

class Genome::Model::Tools::BamQc::SummarizeAsText {
    is => ['Genome::Model::Tools::BamQc::Base'],
    has_input => [
        labels => {
            doc => 'A comma-delimited list of unique labels per directory.',
        },
        directories => {
            doc => 'A comma-delimited list of BamQc directories.',
        },
        output_basename => {
            doc => 'The full path basename to use for output summary files.',
        },
        labels_are_instrument_data_ids => { #this is too long. can't think of anything better.
            is => "Boolean",
            default => 0,
            doc => 'Assume the labels are instrument data ids and query out the library names from the database',
        },
    ],
};

sub _validate_labels_and_directories {
    my $self = shift;

    my @labels = $self->_labels_list;
    my @directories = $self->_directories_list;
    return scalar(@labels) == scalar(@directories);
}

sub _labels_list {
    my $self = shift;

    my $labels_string = $self->labels;
    my @labels = split(',',$labels_string);
    return @labels;
}

sub _directories_list {
    my $self = shift;

    my $directories_string = $self->directories;
    my @directories = split(',',$directories_string);
    return @directories;
}

sub _create_metrics_summary_writer {
    my $self = shift;

    my @summary_headers = qw/
                                LABEL
                                READS
                                READS_ALIGNED
                                PCT_READS_ALIGNED
                                ALIGNED_BASES
                                PCT_CHIMERAS
                                LIBRARY_NAME
                                PCT_DUPLICATION
                                ESTIMATED_LIBRARY_SIZE
                                ERROR_RATE_READ_1
                                ERROR_RATE_READ_2
                                MEDIAN_INSERT_SIZE
                                MEAN_INSERT_SIZE
                                STANDARD_DEVIATION
                                AT_DROPOUT
                                GC_DROPOUT
                            /;
    my $summary_file = $self->output_basename .'-MetricsSummary.tsv';
    if (-e $summary_file) {
        unlink $summary_file;
    }
    my $summary_writer = Genome::Utility::IO::SeparatedValueWriter->create(
        output => $summary_file,
        separator => "\t",
        headers => \@summary_headers,
    );
    return $summary_writer;
}

sub execute {
    my $self = shift;

    $self->_validate_labels_and_directories()
        or die 'Unbalanced input labels and directories!';

    my $summary_writer = $self->_create_metrics_summary_writer();

    my %error_rate_by_position;
    my %is_data;
    my %is_directions;

    my %gc_data;
    my %gc_windows;

    my %qd_data;

    my %qc_data;
    my @labels = $self->_labels_list;
    my @directories = $self->_directories_list;
    for (my $i = 0; $i< scalar(@labels); $i++) {
        my $label = $labels[$i];
        my $output_dir = $directories[$i];

        # Load Alignment Metrics
        my ($as_file) = glob($output_dir .'/*alignment_summary_metrics');
        unless ($as_file) {
            die ('Failed to find Picard alignment_summary_metrics in directory: '. $output_dir);
        }
        my $as_metrics = Genome::Model::Tools::Picard::CollectAlignmentSummaryMetrics->parse_file_into_metrics_hashref($as_file);

        
        # Load MarkDuplicates Metrics
        my ($mrkdup_file) = glob($output_dir .'/*.metrics');
        my $mrkdup_metrics;
        my $lib;
        if ($mrkdup_file) {
            $mrkdup_metrics = Genome::Model::Tools::Picard::MarkDuplicates->parse_file_into_metrics_hashref($mrkdup_file);
            my @libs = sort keys %{$mrkdup_metrics};
            if (scalar(@libs) > 1) {
                die('Unable to handle multiple library BAM files!');
            }
            $lib = $libs[0];
        }
        
        # Load Error Rate Metrics
        my ($error_rate_file) = glob($output_dir .'/*-ErrorRate.tsv');
        my %error_rate_sum;
        if ($error_rate_file) {
            my $error_rate_reader = Genome::Utility::IO::SeparatedValueReader->create(
                input => $error_rate_file,
                separator => "\t",
                ignore_lines_starting_with => '#',
            );
            while (my $error_rate_data = $error_rate_reader->next) {
                if ($error_rate_data->{position} eq 'SUM') {
                    $error_rate_sum{$error_rate_data->{read_end}} = $error_rate_data;
                } else {
                    $error_rate_by_position{$error_rate_data->{read_end}}{$error_rate_data->{position}}{$label} = $error_rate_data->{error_rate};
                }
            }
        }
        
        # Load Insert Size Metrics
        my ($is_file) = glob($output_dir .'/*.insert_size_metrics');
        my $is_metrics = Genome::Model::Tools::Picard::CollectInsertSizeMetrics->parse_file_into_metrics_hashref($is_file);
        
        # Load Insert Size Histogram
        my $is_histo = Genome::Model::Tools::Picard::CollectInsertSizeMetrics->parse_metrics_file_into_histogram_hashref($is_file);

        for my $is_key (keys %{$is_histo}) {
            my $is_size = $is_histo->{$is_key}{insert_size};
            for my $direction (grep {$_ !~ /^insert_size$/} keys %{$is_histo->{$is_key}}) {
                $is_data{$is_size}{$direction}{$label} = $is_histo->{$is_key}{$direction};
                $is_directions{$direction} = 1;
            }
        }
        
        # Load G+C Bias Metrics
        my ($gc_file) = glob($output_dir .'/*-PicardGC_metrics.txt');
        my ($gc_summary) = glob($output_dir .'/*-PicardGC_summary.txt');
        my $gc_metrics = Genome::Model::Tools::Picard::CollectGcBiasMetrics->parse_file_into_metrics_hashref($gc_summary);
        my $gc_data = Genome::Model::Tools::Picard::CollectGcBiasMetrics->parse_file_into_metrics_hashref($gc_file);
        for my $gc_key (keys %{$gc_data}) {
            my $gc_bin = $gc_data->{$gc_key}{GC};
            $gc_data{$gc_bin}{$label}{NORMALIZED_COVERAGE} = $gc_data->{$gc_key}{NORMALIZED_COVERAGE};
            unless ($gc_windows{$gc_bin}) {
                $gc_windows{$gc_bin} = $gc_data->{$gc_key}{WINDOWS};
            } else {
                unless ($gc_windows{$gc_bin} == $gc_data->{$gc_key}{WINDOWS}) {
                    die ($label .' '. $gc_key);
                }
            }
        }

        # Load Quality Distribution
        my ($qd_file) = glob($output_dir .'/*.quality_distribution_metrics');
        my $qd_histo = Genome::Model::Tools::Picard->parse_metrics_file_into_histogram_hashref($qd_file);
        for my $quality_key (keys %{$qd_histo}) {
            my $quality = $qd_histo->{$quality_key}{QUALITY};
            $qd_data{$quality}{$label} = $qd_histo->{$quality_key}{COUNT_OF_Q};
        }

        # Load Quality by Cycle
        my ($qc_file) = glob($output_dir .'/*.quality_by_cycle_metrics');
        my $qc_histo = Genome::Model::Tools::Picard->parse_metrics_file_into_histogram_hashref($qc_file);
        for my $cycle_key (keys %{$qc_histo}) {
            my $cycle = $qc_histo->{$cycle_key}{CYCLE};
            $qc_data{$cycle}{$label} = $qc_histo->{$cycle_key}{MEAN_QUALITY};
        }
        
        # TODO: We need summary metrics per category and/or read direction
        # The below summary metrics only apply to paired-end libraries
        my %summary_data = (
            LABEL => $label,
            READS => $as_metrics->{'CATEGORY-PAIR'}{PF_READS},
            READS_ALIGNED => $as_metrics->{'CATEGORY-PAIR'}{PF_READS_ALIGNED},
            PCT_READS_ALIGNED => $as_metrics->{'CATEGORY-PAIR'}{PCT_PF_READS_ALIGNED},
            ALIGNED_BASES => $as_metrics->{'CATEGORY-PAIR'}{PF_ALIGNED_BASES},
            PCT_CHIMERAS => $as_metrics->{'CATEGORY-PAIR'}{PCT_CHIMERAS},
            ERROR_RATE_READ_1 => $error_rate_sum{1}->{error_rate},
            ERROR_RATE_READ_2 => $error_rate_sum{2}->{error_rate},
            MEDIAN_INSERT_SIZE => $is_metrics->{'PAIR_ORIENTATION-FR'}{MEDIAN_INSERT_SIZE},
            MEAN_INSERT_SIZE => $is_metrics->{'PAIR_ORIENTATION-FR'}{MEAN_INSERT_SIZE},
            STANDARD_DEVIATION => $is_metrics->{'PAIR_ORIENTATION-FR'}{STANDARD_DEVIATION},
            AT_DROPOUT => $gc_metrics->{'WINDOW_SIZE-100'}{AT_DROPOUT},
            GC_DROPOUT => $gc_metrics->{'WINDOW_SIZE-100'}{GC_DROPOUT},
            LIBRARY_NAME => $lib || 'na',
            PCT_DUPLICATION => 'na',
            ESTIMATED_LIBRARY_SIZE => 'na',
        );
        if ($mrkdup_metrics && $lib) {
            $summary_data{PCT_DUPLICATION} = $mrkdup_metrics->{$lib}{PERCENT_DUPLICATION};
            $summary_data{ESTIMATED_LIBRARY_SIZE} = $mrkdup_metrics->{$lib}{ESTIMATED_LIBRARY_SIZE};
        }
        if(!defined($lib) and $self->labels_are_instrument_data_ids) {
            $self->status_message("Looking up the library name");
            my $instrument_data = Genome::InstrumentData->get($label);
            unless($instrument_data) {
                $self->error_message("Unable to find Instrument data for id: $label");
            }
            else {
                my $library_name = $instrument_data->library_name;
                if(defined $library_name) {
                    $summary_data{LIBRARY_NAME} = $library_name;
                }
            }
        }
        $summary_writer->write_one(\%summary_data);
    }
    
    # Error tsv file, wait till the end so we know the maximum position
    for my $read_end (keys %error_rate_by_position) {
        my @error_rate_by_position_headers = ('position',@labels);
        my $error_rate_file = $self->output_basename .'-ErrorRateByPositionRead'. $read_end .'.tsv';
        if (-e $error_rate_file) {
            unlink $error_rate_file;
        }
        my $writer = Genome::Utility::IO::SeparatedValueWriter->create(
            output => $error_rate_file,
            separator => "\t",
            headers => \@error_rate_by_position_headers,
        );
        for my $position (sort {$a <=> $b} keys %{$error_rate_by_position{$read_end}}) {
            my %data = (
                position => $position,
            );
            for my $label (@labels) {
                $data{$label} = $error_rate_by_position{$read_end}{$position}{$label} || 0;
            }
            $writer->write_one(\%data);
        }
    }

    # Write a consolidate histogram of normalized coverage per GC window
    my $gc_summary = $self->output_basename .'-GcBias.tsv';
    if (-e $gc_summary) {
        unlink($gc_summary);
    }
    my @gc_headers = ('GC','WINDOWS',@labels);
    my $gc_data_writer = Genome::Utility::IO::SeparatedValueWriter->create(
        output => $gc_summary,
        separator => "\t",
        headers => \@gc_headers,
    );
    for my $gc_bin (sort {$a <=> $b} keys %gc_windows)  {
        my %data = (
            GC => $gc_bin,
            WINDOWS => $gc_windows{$gc_bin},
        );
        for my $label (@labels) {
            $data{$label} = $gc_data{$gc_bin}{$label}{NORMALIZED_COVERAGE} || 0;
        }
        $gc_data_writer->write_one(\%data);
    }
    
    # Write a consolidate histogram of base quality
    my $qd_summary = $self->output_basename .'-QualityDistribution.tsv';
    if (-e $qd_summary) {
        unlink $qd_summary;
    }
    my @qd_headers = ('QUALITY',@labels);
    my $qd_data_writer = Genome::Utility::IO::SeparatedValueWriter->create(
        output => $qd_summary,
        separator => "\t",
        headers => \@qd_headers,
    ); 
    for my $quality (sort {$a <=> $b} keys %qd_data) {
        my %data = (
            'QUALITY' => $quality
        );
        for my $label (@labels) {
            $data{$label} = $qd_data{$quality}{$label};
        }
        $qd_data_writer->write_one(\%data);
    }

    # Write a consolidate histogram of base quality by cycle
    my $qc_summary = $self->output_basename .'-QualityByCycle.tsv';
    if (-e $qc_summary) {
        unlink $qc_summary;
    }
    my @qc_headers = ('CYCLE',@labels);
    my $qc_data_writer = Genome::Utility::IO::SeparatedValueWriter->create(
        output => $qc_summary,
        separator => "\t",
        headers => \@qc_headers,
    );
    for my $cycle (sort {$a <=> $b} keys %qc_data) {
        my %data = (
            'CYCLE' => $cycle,
        );
        for my $label (@labels) {
            $data{$label} = $qc_data{$cycle}{$label};
        }
        $qc_data_writer->write_one(\%data);
    }

    # Write a consolidate histogram of insert sizes by read orientation/direction
    for my $direction (sort keys %is_directions) {
        my $is_summary = $self->output_basename .'-'. uc($direction) .'-InsertSize.tsv';
        if (-e $is_summary) {
            unlink $is_summary;
        }
        my @is_headers =('INSERT_SIZE',@labels);
        my $is_data_writer = Genome::Utility::IO::SeparatedValueWriter->create(
            output => $is_summary,
            separator => "\t",
            headers => \@is_headers,
        );
        for my $is_size (sort {$a <=> $b} keys %is_data) {
            my %data = (
                INSERT_SIZE  => $is_size,
            );
            for my $label (@labels) {
                $data{$label} = $is_data{$is_size}{$direction}{$label} || 0;
            }
            $is_data_writer->write_one(\%data);
        }
    }
    return 1;
}



1;
