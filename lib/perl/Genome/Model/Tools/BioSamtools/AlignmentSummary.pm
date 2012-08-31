package Genome::Model::Tools::BioSamtools::AlignmentSummary;

use strict;
use warnings;

use Genome;

my %SORT_ORDER = (
    total_bp => 1,
    total_aligned_bp => 2,
    total_unaligned_bp => 3,
    total_duplicate_bp => 4,
    paired_end_bp => 5,
    read_1_bp => 6,
    read_2_bp => 7,
    mapped_paired_end_bp => 8,
    proper_paired_end_bp => 9,
    singleton_bp => 10,
    total_target_aligned_bp => 11,
    unique_target_aligned_bp => 12,
    duplicate_target_aligned_bp => 13,
    total_off_target_aligned_bp => 14,
    unique_off_target_aligned_bp => 15,
    duplicate_off_target_aligned_bp => 16,
);

class Genome::Model::Tools::BioSamtools::AlignmentSummary {
    is => ['Genome::Model::Tools::BioSamtools'],
    has_input => [
        bam_file => {
            is => 'Text',
            doc => 'A BAM format file of alignment data'
        },
        bed_file => {
            is => 'Text',
            doc => 'A BED file of regions of interest to evaluate on/off target alignment',
            is_optional => 1,
        },
        wingspan => {
            is => 'Integer',
            doc => 'A wingspan to add to each region of interest coordinate span',
            default_value => 0,
            is_optional => 1,
        },
        output_directory => {
            is => 'Text',
            doc => 'A directory path used to resolve the output file name.  Required if output_file not defined.  Mostly used for workflow parallelization.',
            is_optional => 1,
        },
    ],
    has_output => [
        output_file => {
            is => 'Text',
            doc => 'A file path to store tab delimited output.  Required if ouput_directory not provided.',
            is_optional => 1,
        },
    ],
    has_param => [
        lsf_queue => {
            doc => 'When run in parallel, the LSF queue to submit jobs to.',
            is_optional => 1,
            default_value => 'apipe',
        },
        lsf_resource => {
            doc => 'When run in parallel, the resource request necessary to run jobs on LSF.',
            is_optional => 1,
            default_value => "-R 'select[type==LINUX64]'",
        },
    ],
};

sub execute {
    my $self = shift;

    my $bam_file = $self->bam_file;
    # resolve the output file but only use it if the param was not defined
    my $resolved_output_file;
    my $output_directory = $self->output_directory;
    if ($output_directory) {
        unless ($output_directory) {
            die('Failed to provide either output_file or output_directory!');
        }
        unless (-d $output_directory) {
            unless (Genome::Sys->create_directory($output_directory)) {
                die('Failed to create output directory: '. $output_directory);
            }
        }
        my ($bam_basename,$bam_dirname,$bam_suffix) = File::Basename::fileparse($bam_file,qw/\.bam/);
        unless(defined($bam_suffix)) {
            die ('Failed to recognize bam_file '. $bam_file .' without bam suffix');
        }
        $resolved_output_file = $output_directory .'/'. $bam_basename;
    } elsif (!defined($self->output_file)) {
        die('Failed to provide either output_file or output_directory!');
    }

    my $bed_file = $self->bed_file;
    my $regions;
    if ($bed_file) {
        my ($bed_basename,$bed_dirname,$bed_suffix) = File::Basename::fileparse($bed_file,qw/\.bed/);
        unless(defined($bed_suffix)) {
            die ('Failed to recognize bed_file '. $bed_file .' without bed suffix');
        }
        $resolved_output_file .= '-'. $bed_basename;
        if (defined($self->wingspan)) {
            $resolved_output_file .= '-wingspan_'. $self->wingspan;
        }
        $regions = Genome::Model::Tools::RefCov::ROI::Bed->create(
            file => $bed_file,
            #Used for fast lookups at the cost of memory, the higher the number the more memory
            region_index_substring => 5,
            wingspan => $self->wingspan,
        );
        unless ($regions) {
            die('Failed to load region file '. $bed_file .'.  Accepted formats are: bed');
        }
    }
    $resolved_output_file .= '-alignment_summary.tsv';
    unless ($self->output_file) {
        $self->output_file($resolved_output_file);
    }
    my $out_fh = Genome::Sys->open_file_for_writing($self->output_file);

    my $refcov_bam  = Genome::Model::Tools::RefCov::Bam->create(bam_file => $bam_file );
    unless ($refcov_bam) {
        die('Failed to load bam file '. $bam_file);
    }
    my $bam  = $refcov_bam->bio_db_bam;
    my $header = $bam->header();

    # Number of reference sequences
    my $targets = $header->n_targets();

    # The reference sequence names in an array ref with indexed positions
    my $target_names = $header->target_name();

    # at the low level API the seq_id/target_name is meaningless
    # cache the target_names in a hash by actual reference sequence name
    # then we can look up the target index on the fly
    my %target_name_index;
    my $i = 0;
    for my $target_name (@{ $target_names }) {
        $target_name_index{$target_name} = $i++;
    }

    # Make sure our index is not off
    unless ($targets == $i) {
        die 'Expected '. $targets .' targets but counted '. $i .' indices';
    }

    my %stats_summary;
    while (my $align = $bam->read1()) {
        my $len = $align->l_qseq;
        $stats_summary{total_bp} += $len;
        my $duplicate;
        my $flag = $align->flag;
        if ($flag & 1024) {
            $duplicate = 1;
            $stats_summary{total_duplicate_bp} += $len;
        }
        unless ($flag & 4) {
            # These methods may be more accurate to calculate the length of the actual alignment
            # $align->cigar2qlen
            # OR
            # my $start = $align->pos + 1;
            # my $end = $align->calend;
            # my $len = $end - $start + 1;
            $stats_summary{total_aligned_bp} += $len;
            if ($regions) {
                my $chrom = $target_names->[$align->tid];
                my $start = $align->pos + 1;
                my $end = $align->calend;
                #  We could get fancy and actually calculate the amount of overlap...
                if ($regions->overlaps_regions($chrom,$start,$end)) {
                    $stats_summary{total_target_aligned_bp} += $len;
                    if ($duplicate) {
                        $stats_summary{duplicate_target_aligned_bp} += $len;
                    } else {
                        $stats_summary{unique_target_aligned_bp} += $len;
                    }
                } else {
                    $stats_summary{total_off_target_aligned_bp} += $len;
                    if ($duplicate) {
                        $stats_summary{duplicate_off_target_aligned_bp} += $len;
                    } else {
                        $stats_summary{unique_off_target_aligned_bp} += $len;
                    }
                }
            }
        } else {
            $stats_summary{total_unaligned_bp} += $len;
        }
        if ($flag & 1) {
            $stats_summary{paired_end_bp} += $len;
            if ($flag & 2) {
                $stats_summary{proper_paired_end_bp} += $len;
                # TODO: it may be useful to calculate the insert size with stdev...
                # Would this only be for proper pairs or all paired-end reads...
                #$align->isize;
            }
            if ($flag & 64)  {
                $stats_summary{read_1_bp} += $len;
            } else {
                $stats_summary{read_2_bp} += $len;
            }
            if (!($flag & 4) && !($flag & 8)) {
                $stats_summary{mapped_paired_end_bp} += $len;
                unless ($flag & 2) {
                    $stats_summary{singleton_bp} += $len;
                }
            }
        }
    }
    my @headers = sort hash_sort_order (keys %stats_summary);
    print $out_fh join("\t", @headers) ."\n";
    print $out_fh join("\t", map { defined $_ ? $_ : '' } map { $stats_summary{$_} } @headers) ."\n";
    return 1;
}

sub hash_sort_order {
    $SORT_ORDER{$a} <=> $SORT_ORDER{$b};
}

1;
