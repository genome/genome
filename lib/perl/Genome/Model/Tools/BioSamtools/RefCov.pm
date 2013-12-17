package Genome::Model::Tools::BioSamtools::RefCov;

use strict;
use warnings;

use Genome;

class Genome::Model::Tools::BioSamtools::RefCov {
    is => ['Genome::Model::Tools::BioSamtools'],
    has_input => [
        bam_file => {
            doc => 'A BAM file, sorted and indexed, containing alignment/read data',
        },
        bed_file => {
            doc => 'The BED format file (tab delimited: chr,start,end,name) file containing annotation or regions of interest.',
        },
        output_directory => {
            doc => 'When run in parallel, this directory will contain all output and intermediate STATS files. Sub-directories will be made for wingspan and min_depth_filter params. Do not define if stats_file is defined.',
            is_optional => 1,
        },
        min_depth_filter => {
            doc => 'The minimum depth at each position to consider coverage.  For more than one, supply a comma delimited list(ie. 1,5,10,15,20)',
            default_value => 1,
            is_optional => 1,
        },
        wingspan => {
            doc => 'A base pair wingspan value to add +/- of the input regions',
            is_optional => 1,
        },
        min_base_quality => {
            doc => 'only consider bases with a minimum phred quality',
            default_value => 0,
            is_optional => 1,
        },
        min_mapping_quality => {
            doc => 'only consider alignments with minimum mapping quality',
            default_value => 0,
            is_optional => 1,
        }
    ],
    has_output => [
        stats_file => {
            doc => 'When run in parallel, do not define.  From the command line this file will contain the output metrics for each region.',
            is_optional => 1,
        },
        final_directory => {
            doc => 'The directory where parallel output is written to',
            is_optional => 1,
        },
    ],
    has_param => [
        lsf_queue => {
            doc => 'When run in parallel, the LSF queue to submit jobs to.',
            is_optional => 1,
            default_value => $ENV{GENOME_LSF_QUEUE_BUILD_WORKER_ALT},
        },
        lsf_resource => {
            doc => 'When run in parallel, the resource request necessary to run jobs on LSF.',
            is_optional => 1,
            default_value => "-R 'select[type==LINUX64]'",
        },
    ],
};

sub help_detail {
'
These commands are setup to run perl v5.10.0 scripts that use Bio-Samtools and require bioperl v1.6.0.  They all require 64-bit architecture.

Output file format(stats_file):
[1] Region Name (column 4 of BED file)
[2] Percent of Reference Bases Covered
[3] Total Number of Reference Bases
[4] Total Number of Covered Bases
[5] Number of Missing Bases
[6] Average Coverage Depth
[7] Standard Deviation Average Coverage Depth
[8] Median Coverage Depth
[9] Number of Gaps
[10] Average Gap Length
[11] Standard Deviation Average Gap Length
[12] Median Gap Length
[13] Min. Depth Filter
[14] Discarded Bases (Min. Depth Filter)
[15] Percent Discarded Bases (Min. Depth Filter)
';
}

sub execute {
    my $self = shift;

    my $output_directory = $self->output_directory;
    my $wingspan = $self->wingspan;
    if ($output_directory) {
        # TODO: Old CoverageStats required running RefCov with multiple Wingspan valude
        # In reality, this create redundancy across the ROI.

        # Moving forward CoverageStatsV2 will create multiple BED files.
        # What do we do with this auto-generated subdirectory....
        if (defined($wingspan)) {
            $output_directory .= '/wingspan_'. $wingspan;
        }
        unless (-d $output_directory){
            unless (Genome::Sys->create_directory($output_directory)) {
                die('Failed to create output directory '. $output_directory);
            }
        }
        $self->final_directory($output_directory);
    }
    unless (defined($self->stats_file)) {
        unless (defined($self->output_directory)) {
            die('Failed to define output_directory or stats_file!');
        }
        my ($bam_basename,$bam_dirname,$bam_suffix) = File::Basename::fileparse($self->bam_file,qw/.bam/);
        unless (defined($bam_suffix)) {
            die('Failed to recognize bam_file '. $self->bam_file .' without bam suffix');
        }
        my ($regions_basename,$regions_dirname,$bed_suffix) = File::Basename::fileparse($self->bed_file,qw/.bed/);
        unless (defined($bed_suffix)) {
            die('Failed to recognize bed_file '. $self->bed_file .' without bed suffix');
        }
        $self->stats_file($self->final_directory .'/'. $bam_basename .'_'. $regions_basename .'_STATS.tsv');
    }

    my $temp_stats_file = Genome::Sys->create_temp_file_path;
    $self->warning_message('Please use \'gmt ref-cov\' instead of \'gmt bio-samtools ref-cov\'!');
    my $cmd = 'genome-perl5.10 -S gmt ref-cov standard --alignment-file-path='. $self->bam_file .' --min-depth-filter='. $self->min_depth_filter .' --roi-file-path='. $self->bed_file .' --stats-file='. $temp_stats_file;
    if (defined $wingspan) {
        $cmd .= ' --wingspan='. $wingspan ;
    }
    if ($self->min_base_quality) {
        $cmd .= ' --min-base-quality='. $self->min_base_quality;
    }
    if ($self->min_mapping_quality) {
        $cmd .= ' --min-mapping-quality='. $self->min_mapping_quality;
    }
    Genome::Sys->shellcmd(
        cmd => $cmd,
        input_files => [$self->bam_file,$self->bed_file],
        output_files => [$temp_stats_file],
    );

    Genome::Sys->copy_file($temp_stats_file, $self->stats_file);

    return 1;
}

1;
