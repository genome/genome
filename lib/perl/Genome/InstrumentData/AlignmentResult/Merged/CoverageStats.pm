package Genome::InstrumentData::AlignmentResult::Merged::CoverageStats;

use strict;
use warnings;

use Genome;
use Sys::Hostname;
use File::Path;

class Genome::InstrumentData::AlignmentResult::Merged::CoverageStats {
    is => ['Genome::SoftwareResult::Stageable'],
    has_input => [
        alignment_result_id => {
            is => 'Number',
            doc => 'ID of the result for the alignment data upon which to run coverage stats',
        },
        region_of_interest_set_id => {
            is => 'Text',
            doc => 'ID of the feature list containing the regions over which to gather coverage stats',
        },
    ],
    has_param => [
        minimum_depths => {
            is => 'Text',
            doc => 'comma-separated list of minimum depths at which to evaluate coverage',
        },
        wingspan_values => {
            is => 'Text',
            doc => 'comma-separated list of wingspan values to add to each region',
        },
        minimum_base_quality => {
            is => 'Text',
            doc => 'A minimum base quality to consider in coverage',
        },
        minimum_mapping_quality => {
            is => 'Text',
            doc => 'A minimum mapping quality to consider in coverage',
        },
        use_short_roi_names => {
            is => 'Boolean',
            doc => 'Whether or not to shorten the names in the BED file for processing',
        },
        roi_track_name => {
            is => 'Text',
            doc => 'For multi-tracked ROI use this named track',
            is_optional => 1,
        },
        merge_contiguous_regions => {
            is => 'Boolean',
            doc => 'Whether or not to merge overlapping/adjoining regions before analysis',
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
            is => 'Genome::InstrumentData::AlignedBamResult',
            id_by => 'alignment_result_id',
            doc => 'the alignment data upon which to run coverage stats',
        },
        region_of_interest_set => {
            is => 'Genome::FeatureList',
            id_by => 'region_of_interest_set_id',
            doc => 'regions over which to gather coverage stats',
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
    my $base_dir = sprintf("coverage-%s-%s-%s-%s", $hostname, $user, $$, $self->id);

    # TODO: the first subdir is actually specified by the disk management system.
    my $directory = join('/', 'build_merged_alignments','coverage_stats',$base_dir);
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
    return 'coverage-stats';
}

sub _prepare_staging_directory {
    my $self = shift;

    return $self->temp_staging_directory if ($self->temp_staging_directory);

    unless($self->output_dir) {
        $self->_prepare_output_directory;
    }

    #Stage to network disk because of inner workflow
    my $staging_tempdir = File::Temp->newdir(
        $self->_working_dir_prefix . '-staging-XXXXX',
        DIR     => $self->output_dir,
        CLEANUP => 1,
    );

    $self->temp_staging_directory($staging_tempdir);
}

sub create {
    my $class = shift;
    my $self = $class->SUPER::create(@_);
    return unless ($self);

    $self->_prepare_staging_directory;

    my $bed_file = $self->_dump_bed_file;
    my $bam_file = $self->alignment_result->bam_file;

    die $self->error_message("Bed File ($bed_file) is missing") unless -s $bed_file;
    die $self->error_message("Bam File ($bam_file) is missing") unless -s $bam_file;

    my $log_dir = $self->log_directory;
    unless($log_dir) {
        $log_dir = '' . $self->temp_staging_directory;
    }
    $self->_log_directory($log_dir);

    my %coverage_params = (
        output_directory => '' . $self->temp_staging_directory,
        log_directory => $log_dir,
        bed_file => $bed_file,
        bam_file => $bam_file,
        minimum_depths => $self->minimum_depths,
        wingspan_values => $self->wingspan_values,
        minimum_base_quality => $self->minimum_base_quality,
        minimum_mapping_quality => $self->minimum_mapping_quality,
    );

    unless($] > 5.010) {
        #need to shell out to a newer perl #TODO remove this once 5.10 transition complete
        my $cmd = 'genome-perl5.10 -S gmt bio-samtools coverage-stats ';
        while (my ($key, $value) = (each %coverage_params)) {
            $key =~ s/_/-/g;
            $cmd .= " --$key=$value";
        }

        Genome::Sys->shellcmd(
            cmd => $cmd,
            input_files => [$bed_file, $bam_file],
        );
    } else {
        my $cmd = Genome::Model::Tools::BioSamtools::CoverageStats->create(%coverage_params);
        unless($cmd->execute) {
            die('Failed to run coverage stats tool');
        }
    }

    $self->_promote_data;
    $self->_reallocate_disk_allocation;

    $self->_generate_metrics;

    #quick sanity check--all wingspans should run on the same number of total reads
    my @wingspan_values = split(',', $self->wingspan_values);
    if(@wingspan_values > 1) {
        my $total_bp = 0;

        my $alignment_summary = $self->alignment_summary_hash_ref;
        for my $wingspan (@wingspan_values) {
            $total_bp ||= $alignment_summary->{$wingspan}{total_bp}; #initialize first time through loop

            unless($total_bp eq $alignment_summary->{$wingspan}{total_bp}) {
                die($self->error_message('Inconsistent results produced.'));
            }
        }
    }


    return $self;
}

sub _dump_bed_file {
    my $self = shift;

    my $roi_set = $self->region_of_interest_set;
    return unless $roi_set;

    my $alt_reference;
    my $reference = $self->alignment_result->reference_build;
    unless($reference->is_compatible_with($roi_set->reference)) {
        $alt_reference = $reference;
    }
    my $merge_status = $self->merge_contiguous_regions;
    my $use_short_names = $self->use_short_roi_names;

    my $bed_file_path = $self->temp_staging_directory .'/'. $roi_set->id .'.bed';
    unless (-e $bed_file_path) {
        my %dump_params = (
            feature_list => $roi_set,
            output_path => $bed_file_path,
            alternate_reference => $alt_reference,
            merge => $merge_status,
            short_name => $use_short_names,
        );
        if ($self->roi_track_name) {
            $dump_params{track_name} = $self->roi_track_name;
        }
        my $dump_command = Genome::FeatureList::Command::DumpMergedList->create(%dump_params);
        unless ($dump_command->execute) {
            die('Failed to print bed file to path '. $bed_file_path);
        }
    }

    return $bed_file_path;
}

sub _generate_metrics {
    my $self = shift;

    $self->_generate_alignment_summary_metrics('');;
    $self->_generate_alignment_summary_metrics('-v2');
    $self->_generate_coverage_stats_summary_metrics;

    return 1;
}

sub alignment_summary_file {
    my ($self,$wingspan,$version_suffix) = @_;
    unless (defined($wingspan)) {
        die('Must provide wingspan_value to method alignment_summary_file in '. __PACKAGE__);
    }

    # $version suffix should either be an empty string or something starting with '-'
    if ($version_suffix) {
        if (substr($version_suffix, 0, 1) ne '-') {
            $version_suffix = '-' . $version_suffix;
        }
    } else {
        $version_suffix = '';
    }

    my @as_files = glob(sprintf('%s/*wingspan_%s-alignment_summary%s.tsv',
                                $self->output_dir,
                                $wingspan,
                                $version_suffix));
    unless (@as_files) {
        return;
    }
    unless (scalar(@as_files) == 1) {
        die("Found multiple stats files:\n". join("\n",@as_files));
    }
    return $as_files[0];
}

sub _generate_alignment_summary_metrics {
    my $self = shift;
    my $version_suffix = shift;

    my $version ||= '';

    my @wingspan_values = split(',', $self->wingspan_values);

    for my $wingspan (@wingspan_values) {
        my $data;

        my $as_file = $self->alignment_summary_file($wingspan, $version_suffix);
        my $reader = Genome::Utility::IO::SeparatedValueReader->create(
            separator => "\t",
            input => $as_file,
        );
        unless ($reader) {
            $self->error_message('Can not create SeparatedValueReader for input file '. $as_file);
            die($self->error_message);
        }
        $data = $reader->next;
        $reader->input->close;
        # Calculate percentages

        # percent aligned
        $data->{percent_aligned} = $data->{total_bp}
                                    ? sprintf("%.02f",(($data->{total_aligned_bp} / $data->{total_bp}) * 100))
                                    : 0;

        # duplication rate
        $data->{percent_duplicates} = $data->{total_aligned_bp}
                                    ? sprintf("%.03f",(($data->{total_duplicate_bp} / $data->{total_aligned_bp}) * 100))
                                    : 0;

        # on-target alignment
        $data->{percent_target_aligned} = $data->{total_aligned_bp}
                                    ? sprintf("%.02f",(($data->{total_target_aligned_bp} / $data->{total_aligned_bp}) * 100))
                                    : 0;

        # on-target duplicates
        $data->{percent_target_duplicates} = $data->{total_target_aligned_bp}
                                    ? sprintf("%.02f",(($data->{duplicate_target_aligned_bp} / $data->{total_target_aligned_bp}) * 100))
                                    : 0;

        # off-target alignment
        $data->{percent_off_target_aligned} = $data->{total_aligned_bp}
                                    ? sprintf("%.02f",(($data->{total_off_target_aligned_bp} / $data->{total_aligned_bp}) * 100))
                                    : 0;

        # off-target duplicates
        $data->{percent_off_target_duplicates} = $data->{total_off_target_aligned_bp}
                                    ? sprintf("%.02f",(($data->{duplicate_off_target_aligned_bp} / $data->{total_off_target_aligned_bp}) * 100))
                                    : 0;

        for my $key (keys %$data) {
            my $metric_key = sprintf('alignment%s-wingspan_%s_%s',
                                    $version_suffix, $wingspan, $key);
            $self->add_metric(metric_name => $metric_key, metric_value => $data->{$key});
        }
    }

    return 1;
}

sub alignment_summary_hash_ref {
    my $self = shift;

    unless ($self->{_alignment_summary_hash_ref}) {
        my @wingspan_values = split(',', $self->wingspan_values);
        my %alignment_summary;
        for my $wingspan (@wingspan_values) {
            my $alignment_key_basename = 'alignment-wingspan_'. $wingspan;
            my @all_metrics = $self->metrics;
            my @metrics = grep { $_->metric_name =~ /^$alignment_key_basename/ } @all_metrics;
            my $data;

            if (@metrics) {
                for my $metric (@metrics) {
                    my $metric_name = $metric->metric_name;
                    my $alignment_key_regex = $alignment_key_basename .'_(\S+)';
                    unless ($metric_name =~ /^$alignment_key_regex/) {
                        die('Failed to parse alignment metric name '. $metric_name);
                    }
                    my $key = $1;
                    $data->{$key} = $metric->metric_value;
                }
            } else {
                die('No metrics found for this result!');
            }
            $alignment_summary{$wingspan} = $data;
        }
        $self->{_alignment_summary_hash_ref} = \%alignment_summary;
    }
    return $self->{_alignment_summary_hash_ref};
}


sub alignment_summary_v2_hash_ref {
    my $self = shift;

    unless ($self->{_alignment_summary_v2_hash_ref}) {
        my @wingspan_values = split(',', $self->wingspan_values);
        my %alignment_summary;
        for my $wingspan (@wingspan_values) {
            my $alignment_key_basename = 'alignment-v2-wingspan_'. $wingspan;
            my @all_metrics = $self->metrics;
            my @metrics = grep { $_->metric_name =~ /^$alignment_key_basename/ } @all_metrics;
            my $data;

            if (@metrics) {
                for my $metric (@metrics) {
                    my $metric_name = $metric->metric_name;
                    my $alignment_key_regex = $alignment_key_basename .'_(\S+)';
                    unless ($metric_name =~ /^$alignment_key_regex/) {
                        die('Failed to parse alignment metric name '. $metric_name);
                    }
                    my $key = $1;
                    $data->{$key} = $metric->metric_value;
                }
            } else {
                #die('No metrics found for this result!');
                return {};
            }
            $alignment_summary{$wingspan} = $data;
        }
        $self->{_alignment_summary_v2_hash_ref} = \%alignment_summary;
    }
    return $self->{_alignment_summary_v2_hash_ref};
}

sub coverage_stats_directory_path {
    my ($self,$wingspan) = @_;
    return $self->output_dir .'/wingspan_'. $wingspan;
}

sub stats_file {
    my $self = shift;
    my $wingspan = shift;
    unless (defined($wingspan)) {
        die('Must provide wingspan_value to method coverage_stats_file in '. __PACKAGE__);
    }
    my $coverage_stats_directory = $self->coverage_stats_directory_path($wingspan);
    my @stats_files = glob($coverage_stats_directory.'/*STATS.tsv');
    unless (@stats_files) {
        return;
    }
    unless (scalar(@stats_files) == 1) {
        die("Found multiple stats files:\n". join("\n",@stats_files));
    }
    return $stats_files[0];
}

sub coverage_stats_hash_ref {
    my $self = shift;
    unless ($self->{_coverage_stats_hash_ref}) {
        my @headers = qw/name pc_covered length covered_bp uncovered_bp mean_depth stdev_mean_depth median_depth gaps mean_gap_length stdev_gap_length median_gap_length minimum_depth minimum_depth_discarded_bp pc_minimum_depth_discarded_bp/;

        my %stats;
        my @wingspan_values = split(',', $self->wingspan_values);
        for my $wingspan (@wingspan_values) {
            my $stats_file = $self->stats_file($wingspan);
            my $reader = Genome::Utility::IO::SeparatedValueReader->create(
                separator => "\t",
                input => $stats_file,
                #TODO: Add headers to the stats file
                headers => \@headers,
            );
            unless ($reader) {
                $self->error_message('Can not create SeparatedValueReader for file '. $stats_file);
                die $self->error_message;
            }
            while (my $data = $reader->next) {
                push @{$stats{$wingspan}{$data->{name}}}, $data;
            }
            $reader->input->close;
        }
        $self->{_coverage_stats_hash_ref} = \%stats;
    }
    return $self->{_coverage_stats_hash_ref};
}

sub coverage_stats_summary_file {
    my ($self,$wingspan) = @_;
    unless (defined($wingspan)) {
        die('Must provide wingspan_value to method coverage_stats_file in '. __PACKAGE__);
    }
    my $glob_string = $self->coverage_stats_directory_path($wingspan) .'/*STATS.txt';
    my @stats_files = glob($glob_string);
    unless (@stats_files) {
        $self->error_message('Failed to find coverage stats summary file like '. $glob_string);
        die($self->error_message);
    }
    unless (scalar(@stats_files) == 1) {
        die("Found multiple stats summary files:\n". join("\n",@stats_files));
    }
    return $stats_files[0];
}

sub _generate_coverage_stats_summary_metrics {
    my $self = shift;
    my @wingspan_values = split(',', $self->wingspan_values);

    for my $wingspan (@wingspan_values) {
        my $stats_summary = $self->coverage_stats_summary_file($wingspan);
        unless ($stats_summary) {
            $self->error_message('Failed to find coverage stats summary file for wingspan '. $wingspan);
            die($self->error_message);
        }
        my $reader = Genome::Utility::IO::SeparatedValueReader->create(
            separator => "\t",
            input => $stats_summary,
        );
        unless ($reader) {
            $self->error_message('Can not create SeparatedValueReader for file '. $stats_summary);
            die $self->error_message;
        }
        while (my $data = $reader->next) {
            # record stats as build metrics
            for my $key (keys %$data) {
                my $metric_key = join('_', 'coverage-wingspan', $wingspan, $data->{'minimum_depth'}, $key);
                $self->add_metric(metric_name => $metric_key, metric_value => $data->{$key});
            }
        }
        $reader->input->close;
    }

    return 1;
}

sub coverage_stats_summary_hash_ref {
    my $self = shift;
    unless ($self->{_coverage_stats_summary_hash_ref}) {
        my %stats_summary;
        my @wingspan_values = split(',', $self->wingspan_values);
        for my $wingspan (@wingspan_values) {
            my $key_basename = 'coverage-wingspan_'. $wingspan;
            my $min_depth_key_regex = $key_basename .'_(\d+)';
            my @all_metrics = $self->metrics;
            my @metrics = grep { $_->metric_name =~ /^$min_depth_key_regex/ } @all_metrics;

            if (@metrics) {
                for my $metric (@metrics) {
                    my $metric_name = $metric->metric_name;
                    my $coverage_key_regex = $min_depth_key_regex .'_(\S+)';
                    unless ($metric_name =~ /^$coverage_key_regex/) {
                        die('Failed to parse alignment metric name '. $metric_name);
                    }
                    my $min_depth = $1;
                    my $key = $2;
                    $stats_summary{$wingspan}{$min_depth}->{$key} = $metric->metric_value;
                }
            } else {
                die('No metrics found for this result!');
            }
        }
        $self->{_coverage_stats_summary_hash_ref} = \%stats_summary;
    }
    return $self->{_coverage_stats_summary_hash_ref};
}


1;
