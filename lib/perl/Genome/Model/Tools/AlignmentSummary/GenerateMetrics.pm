package Genome::Model::Tools::AlignmentSummary::GenerateMetrics;

use strict;
use warnings;

use Genome;

class Genome::Model::Tools::AlignmentSummary::GenerateMetrics {
    is => 'Genome::Model::Tools::AlignmentSummary::Base',
    has_input => [
        alignment_file_path => {
            doc => 'The aligned file in SAM/BAM format.',
        },
        alignment_file_format => {
            is => 'Text',
            doc => 'The format of the alignment file.',
            valid_values => ['sam','bam'],
            default_value => 'bam',
            is_optional => 1,
        },
        roi_file_path => {
            doc => 'The sorted file with regions of interest.',
            is_optional => 1,
        },
        roi_file_format => {
            is => 'Text',
            doc => 'The format of the sorte annotation file.',
            valid_values => ['bed'],
            default_value => 'bed',
            is_optional => 1,
        },
        verbosity => {
            is => 'Text',
            doc => 'The amount of stdout and stderr to generate',
            valid_values => ['full','no_progress_bars','quiet'],
            default_value => 'quiet',
            is_optional => 1,
        },
        ignore_cigar_md_errors => {
            is => 'Boolean',
            doc => 'rather than exiting upon encountering a CIGAR / MD string error, abort per-bp statistic collection for the current read and move to the next read',
            default_value => 0,
            is_optional => 1,
        },
        no_strict_bases => {
            is => 'Boolean',
            doc => 'accept a-Z as valid bases',
            default_value => 0,
            is_optional => 1,
        },
        output_file => {
            is => 'Text',
            doc => 'The output file',
            is_output => 1,
            is_optional => 1,
        },
        output_directory => {
            is => 'Text',
            doc => 'The output directory when run in parallel via Workflow',
            is_optional => 1,
        },
    ],
};

sub execute {
    my $self = shift;

    my $cmd = $self->path_for_version($self->use_version) .' --'. $self->alignment_file_format .' '. $self->alignment_file_path;
    my @input_files = ($self->alignment_file_path);

    if ($self->verbosity eq 'quiet') {
        $cmd .= ' --quiet';
    } elsif ($self->verbosity eq 'no_progress_bars') {
        $cmd .= ' --no-progress-bars';
    } elsif ($self->verbosity eq 'full') {
        # DO NOTHING
    } else {
        die('Unknown value for verbosity : '. $self->verbosity);
    }
    
    if ($self->roi_file_path) {
        push @input_files, $self->roi_file_path;
        $cmd .= ' --'. $self->roi_file_format .' '. $self->roi_file_path;
    }
    if ($self->ignore_cigar_md_errors) {
        $cmd .= ' --ignore-cigar-md-errors';
    }
    if ($self->no_strict_bases) {
        $cmd .= ' --no-strict-bases';
    }
    
    my $resolved_output_file;
    my $output_directory = $self->output_directory;
    if ($output_directory) {
        unless (-d $output_directory) {
            unless (Genome::Sys->create_directory($output_directory)) {
                die('Failed to create output directory: '. $output_directory);
            }
        }
        my $suffix_re = '.'. $self->alignment_file_format;
        my ($basename,$dirname,$suffix) = File::Basename::fileparse($self->alignment_file_path,$suffix_re);
        unless(defined($suffix)) {
            die ('Failed to recognize alignment file '. $self->alignment_file_path .' without '. $suffix_re .' suffix');
        }
        $resolved_output_file = $output_directory .'/'. $basename;
        my $bed_file = $self->roi_file_path;
        my ($bed_basename,$bed_dirname,$bed_suffix) = File::Basename::fileparse($bed_file,qw/\.bed/);
        unless(defined($bed_suffix)) {
            die ('Failed to recognize bed_file '. $bed_file .' without bed suffix');
        }
        $resolved_output_file .= '-'. $bed_basename .'-alignment_summary.tsv';
        $self->output_file($resolved_output_file);
    }
    unless ($self->output_file) {
        die('Failed to provide output_file or output_directory!');
    }
    $cmd .= ' > '. $self->output_file;
    my @output_files = ($self->output_file);

    $self->run_in_bash_env(
        cmd => $cmd,
        input_files => \@input_files,
        output_files => \@output_files,
    );

    return 1;
}
