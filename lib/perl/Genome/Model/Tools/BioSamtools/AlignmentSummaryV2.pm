package Genome::Model::Tools::BioSamtools::AlignmentSummaryV2;

use strict;
use warnings;

use Genome;

class Genome::Model::Tools::BioSamtools::AlignmentSummaryV2 {
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
};

sub execute {
    my $self = shift;

    eval "require Bio::DB::Sam";
    if ($@) {
        die "Can't use Bio::DB::Sam: $@";
    }
    import Bio::DB::Sam;
    eval "require Bio::DB::Sam::Constants";
    if ($@) {
        die "Can't use Bio::DB::Sam::Constants: $@";
    }
    import Bio::DB::Sam::Constants;

    my $output_file = $self->_resolve_output_file();

    my $bam_reader = $self->_open_sorted_bam();

    $self->_verify_target_index($bam_reader);

    my $temp_bed_directory = $self->_prepare_split_sorted_bed_files($self->wingspan);

    my $bed_reader = $self->_create_bed_reader($bam_reader, $temp_bed_directory);

    $self->debug_message("Summarizing alignments");
    my $stats_summary = $self->summarize_alignments($bam_reader, $bed_reader);

    $self->_save_summary($stats_summary, $output_file);

    return 1;
}

sub _create_bed_reader {
    my($self, $bam_reader, $sorted_bed_directory) = @_;

    my $target_names = $bam_reader->header->target_name();  # converts alignment tid to chrom name
    my $bed_fh;                 # Filehandle for the open bed file
    my $working_tid = \$bed_fh; # tid/chromosome of the open region file
    my($next_range_start, $next_range_stop);

    my $open_range_file = sub {
        my $new_tid = shift;

        $working_tid = $new_tid;
        my $bed_file = $self->_sorted_bed_file_name_for_chrom($sorted_bed_directory, $target_names->[$new_tid]);
        unless (-f $bed_file) {
            # If the file doesn't exist, then there's no regions for it, meaning
            # this alignment has no intersections
            return;
        }
        $bed_fh = IO::File->new($bed_file, 'r');
        $bed_fh || die "Can't open sorted bed file for chromosome " . $target_names->[$new_tid] . " for reading: $!";
    };

    my $read_next_range_from_file = sub {
        if ($_[0] != $working_tid) {
            $open_range_file->($_[0]);
                $working_tid = $_[0];
        $next_range_start = undef;
            $next_range_stop = undef;
        }
        return unless ($bed_fh);   # file handle closed after a previous end-of-file

        my $line = <$bed_fh>;
        if (defined($line) and length($line)) {
            (undef,$next_range_start, $next_range_stop) = split(' ',$line);
        } else {
            # must be end-of-file
            $bed_fh = undef;
            $next_range_start = undef;
        }
        return ($next_range_start, $next_range_stop);
    };

    return $read_next_range_from_file;
}


sub summarize_alignments {
    my($self, $bam_reader, $bed_reader) = @_;

    #$DB::single=1;

    # The statistics we're collecting
    my $total_bp = 0;
    my $total_duplicate_bp = 0;
    my $total_unaligned_bp = 0;
    my $total_aligned_bp = 0;
    my $total_target_aligned_bp = 0;
    my $total_off_target_aligned_bp = 0;
    my $total_soft_clipped_bp = 0;
    my $total_hard_clipped_bp = 0;
    my $total_skipped_bp = 0;
    my $total_deleted_bp = 0;
    my $total_inserted_bp = 0;
    my $duplicate_target_aligned_bp = 0;
    my $duplicate_off_target_aligned_bp = 0;
    my $unique_target_aligned_bp = 0;
    my $unique_off_target_aligned_bp;
    my $paired_end_bp = 0;
    my $proper_paired_end_bp = 0;
    my $read_1_bp = 0;
    my $read_2_bp = 0;
    my $mapped_paired_end_bp = 0;
    my $singleton_bp = 0;

    my $loaded_ranges = []; # list of ranges in the vicinity of the last alignment

    # This closure returns true if this cigar part is eligible to intersect an roi (M or I)
    # It's only used inside the main loop, but is located here to avoid the time required to
    # reassign the same coderef over and over
    my $next_cigar;
    my $cigar_start;   # holds the start/stop for the current cigar segment.  As we iterate
    my $cigar_stop;    # through the segments, we increment these based on past segment lengths
    my $cigar_len;     # length of this particular cigar segment
    my $align_stop;

    my $preprocess_cigar_part = sub {

        $cigar_start = $cigar_stop;
        my $cigar_op  = $next_cigar & &Bio::DB::Sam::Constants::BAM_CIGAR_MASK;
        $cigar_len = $next_cigar >> &Bio::DB::Sam::Constants::BAM_CIGAR_SHIFT;

        if ($cigar_op == 4) {
            # soft clipped: S
            $total_bp += $cigar_len;
            $total_soft_clipped_bp  += $cigar_len;
            $total_unaligned_bp     += $cigar_len;
            return;

        } elsif ($cigar_op == 3) {
            # Skipped: N
            #$total_bp += $cigar_len;
            $total_skipped_bp       += $cigar_len;
            #$total_unaligned_bp     += $cigar_len;
            return

        } elsif ($cigar_op == 5) {
            # Hard clipped: H
            #$total_bp += $cigar_len;
            $total_hard_clipped_bp  += $cigar_len;
            #$total_unaligned_bp     += $cigar_len;
            return

        } elsif ($cigar_op == 2) {
            # Deleted: D
            #$total_bp += $cigar_len;
            $total_deleted_bp  += $cigar_len;
            #$total_unaligned_bp     += $cigar_len;
            return

        } elsif ($cigar_op == 1) {
            # Inserted: I
            $total_inserted_bp += $cigar_len;
            $align_stop += $cigar_len;  # Inserts need to push the alignment stop further down
            #$total_bp += $cigar_len;
            #$total_aligned_bp      += $cigar_len;

        } elsif ($cigar_op != 0) {
            # not a match/mismatch: M
            # Something else?
            die "Don't know what to do with cigar op $cigar_op";
            return;
        }
        $cigar_stop = $cigar_start + $cigar_len;
        return 1;
    };

    my ($next_range_start,$next_range_stop);
    # The main loop

    ALIGNMENT:
    while (my $align = $bam_reader->read1()) {
        #$DB::single=1 if ($align->tid == 0 and $align->pos == 1000 and $align->calend == 1100);

        my $flag = $align->flag;
        my $is_duplicate = $flag & 1024;

        my $alignment_len = $align->l_qseq;

        if ($flag & 4) {
            # entire read is unaligned
            $total_unaligned_bp += $alignment_len;
            $total_bp += $alignment_len;

        } else {
            # it was aligned.  See if it intersects
            my $align_tid = $align->tid;
            my $align_start = $align->pos;
            #my $align_stop = $align->calend + 1;
            $align_stop = $align->calend + 1;

            unless (defined $next_range_start) {
                ($next_range_start, $next_range_stop) = $bed_reader->($align_tid);
            }

            my @intersections;

            # Keep any ranges that cross or are after the alignment
            my @keep_ranges;
            foreach my $range ( @$loaded_ranges ) {
                if ($align_start < $range->[1]) {
                    push @keep_ranges, $range;

                    if ($align_stop >= $range->[0]) {
                        push @intersections, $range;
                    }
                }
            }

            # read in new ranges
            while (defined($next_range_start) and $next_range_start <= $align_stop) {
                my $next_range = [$next_range_start, $next_range_stop];
                push @keep_ranges, $next_range;

                if ($align_start < $next_range_stop) {
                    push @intersections, $next_range;
                }
                ($next_range_start, $next_range_stop) = $bed_reader->($align_tid);
            }

            $loaded_ranges = \@keep_ranges;

            $cigar_start = $align_start;
            $cigar_stop = $align_start;
            if (@intersections) {
                # @intersections now contains a list of ranges that overlap this entire alignment
                #
                # We also need to parse the cigar string of the alignment into smaller segments
                # and consider how each of those pieces intersect the ranges.  Different cigar segment
                # types are counted in different ways.
                #
                # Count up how much overlap there is between the alignment parts we are keeping and the
                # intersecing regions


                my $cigar_parts = $align->cigar;


                my $next_cigar_idx = 0;
                my $next_roi;
                my $next_roi_idx = 0;
                my $current_position = $cigar_start;

                WHILE_IN_READ:
                while($current_position < $align_stop) {
                    if (!$next_cigar or ($current_position >= $cigar_stop)) {
                        $next_cigar = $cigar_parts->[ $next_cigar_idx++ ];
                        last WHILE_IN_READ unless $next_cigar;
                        unless (&$preprocess_cigar_part) {
                            $next_cigar = undef;
                            next WHILE_IN_READ;
                        }
                        $total_bp += $cigar_len;
                        # This counts the entire cigar segment as aligned.
                        # The C++ tool looks at the match-descriptior string to make the
                        # distinction between matching the reference (aligned) and not-matching (unaligned)
                        $total_aligned_bp      += $cigar_len;
                        if ($is_duplicate) {
                            $total_duplicate_bp += $cigar_len;
                        }
                    }
                    if (!$next_roi or ($current_position >= $next_roi->[1])) {
                        $next_roi   = $intersections[ $next_roi_idx++ ];
                    }

                    my $is_in_roi = ($next_roi and $current_position >= $next_roi->[0] and $current_position < $next_roi->[1]);
                    my $is_in_read = ($next_cigar and $current_position >= $cigar_start and $current_position < $cigar_stop);

                    # Figure out how far to advance
                    my $next_position;
                    if ($is_in_read) {
                        if ($is_in_roi) {
                            $next_position = ($next_roi->[1] < $cigar_stop) ? $next_roi->[1] : $cigar_stop;
                        } elsif ($next_roi) {
                            $next_position = ($next_roi->[0] < $cigar_stop) ? $next_roi->[0] : $cigar_stop;
                        } else {
                            $next_position = $cigar_stop;
                        }
                    } elsif ($is_in_roi) {
                        # In the roi, but not in the read, advance to the next position
                        $current_position = ($next_roi->[1] < $cigar_start) ? $next_roi->[1] : $cigar_start;
                        next WHILE_IN_READ;
                    } else {
                        $current_position = $cigar_start;
                        next WHILE_IN_READ;
                    }

                    my $covered_length = $next_position - $current_position;
                    $current_position = $next_position;

                    if ($is_in_roi) {
                        if ($is_in_read) {
                            $total_target_aligned_bp += $covered_length;
                            if ($is_duplicate) {
                                $duplicate_target_aligned_bp += $covered_length;
                            } else {
                                $unique_target_aligned_bp += $covered_length;
                            }
                        }
                    } elsif ($is_in_read) {
                        $total_off_target_aligned_bp += $covered_length;
                        if ($is_duplicate) {
                            $duplicate_off_target_aligned_bp += $covered_length;
                        } else {
                            $unique_off_target_aligned_bp += $covered_length;
                        }
                    }
                } # end while $current_position < $align_end

            } else {
                # The raw alignment did not intersect any regions
                my $cigar_parts = $align->cigar;
                #foreach $next_cigar ( @$cigar_parts ) {
                my $cigar_parts_count = @$cigar_parts;
                for (my $i = 0; $i < $cigar_parts_count; $i++) {
                    $next_cigar = $cigar_parts->[$i];
#$DB::single=1 unless $next_cigar;
                    if (&$preprocess_cigar_part) {
                        $total_bp += $cigar_len;
                        $total_aligned_bp += $cigar_len;
                        $total_off_target_aligned_bp += $cigar_len;
                        if ($is_duplicate) {
                            $total_duplicate_bp += $cigar_len;
                            $duplicate_off_target_aligned_bp += $cigar_len;
                        } else {
                            $unique_off_target_aligned_bp += $cigar_len;
                        }
                    }
                }
            }

        } # end read was aligned

        if ($flag & 1) {
            $paired_end_bp += $alignment_len;

            if ($flag & 2) {
                $proper_paired_end_bp += $alignment_len;
            }

            if ($flag & 64) {
                $read_1_bp += $alignment_len;
            } else {
                $read_2_bp += $alignment_len;
            }

            #if (!($flag & 4) && !($flag & 8)) {
            #unless (($flag & 4) || ($flag & 8)) {
            unless ($flag & 12) {
                $mapped_paired_end_bp += $alignment_len;
                unless ($flag & 2) {
                    $singleton_bp += $alignment_len;
                }
            }
        }

        # These checks are to make sure all the numbers work out consistently.
        # If you make changes to the algorithm, enable these checks during testing
        # and then disable them again when it's working
        #if (($total_target_aligned_bp + $total_off_target_aligned_bp) != $total_aligned_bp) {
        #    printf("Alignment tid %d start %d end %d makes the on/off target counters off\n", $align->tid, $align->pos, $align->calend);
        #    print "difference is ", $total_aligned_bp - $total_target_aligned_bp - $total_off_target_aligned_bp,"\n";
        #    exit;
        #}
        #if (($total_aligned_bp + $total_unaligned_bp) != $total_bp) {
        #    printf("Alignment tid %d start %d end %d makes the total aligned/unaligned counters off\n", $align->tid, $align->pos, $align->calend);
        #    print "difference is ", $total_bp - $total_aligned_bp - $total_unaligned_bp,"\n";
        #    exit;
        #}
        #if (($duplicate_target_aligned_bp + $duplicate_off_target_aligned_bp) != $total_duplicate_bp) {
        #    printf("Alignment tid %d start %d end %d makes the on/off duplicate target counters off\n", $align->tid, $align->pos, $align->calend);
        #    print "difference is ", $total_duplicate_bp - $duplicate_target_aligned_bp - $duplicate_off_target_aligned_bp,"\n";
        #    exit;
        #}


    } # end while next aligment

    return {
        'total_bp'                        => $total_bp,
        'total_duplicate_bp'              => $total_duplicate_bp,
        'total_unaligned_bp'              => $total_unaligned_bp,
        'total_aligned_bp'                => $total_aligned_bp,
        'total_target_aligned_bp'         => $total_target_aligned_bp,
        'total_off_target_aligned_bp'     => $total_off_target_aligned_bp,
        'total_soft_clipped_bp'           => $total_soft_clipped_bp,
        'duplicate_target_aligned_bp'     => $duplicate_target_aligned_bp,
        'duplicate_off_target_aligned_bp' => $duplicate_off_target_aligned_bp,
        'unique_target_aligned_bp'        => $unique_target_aligned_bp,
        'unique_off_target_aligned_bp'    => $unique_off_target_aligned_bp,
        'paired_end_bp'                   => $paired_end_bp,
        'proper_paired_end_bp'            => $proper_paired_end_bp,
        'read_1_bp'                       => $read_1_bp,
        'read_2_bp'                       => $read_2_bp,
        'mapped_paired_end_bp'            => $mapped_paired_end_bp,
        'singleton_bp'                    => $singleton_bp,
    };
}

 
sub _verify_target_index {
    my($self, $bam_reader) = @_;

    $self->debug_message('Verifying target index');

    my $header = $bam_reader->header();

    # Number of reference sequences
    my $targets = $header->n_targets();

    # The reference sequence names in an array ref with indexed positions
    my $target_names = $header->target_name();

    unless( $targets == @{ $target_names }) {
        die "Expected $targets targets, but counted ". scalar(@$target_names). "indices";
    }
    return 1;
}


        
sub _open_sorted_bam {
    my $self = shift;

    if (-f $self->bam_file && -r _) {
        my $bam_reader = Bio::DB::Bam->open( $self->bam_file );
        my $bam_header = $bam_reader->header;
        if ($bam_header->text =~ m/SO:coordinate/i) {
            $self->debug_message("The bam file is already sorted");
            return $bam_reader; 
        }
    } else {
        die "Can't open bam file ".$self->bam_file.": file does not exist or is not readable";
    }

    $self->debug_message('Sorting input bam file');

    my($sort_fh, $sorted_bam) = File::Temp::tempfile('alignment_summary_sorted_XXXXXX', UNLINK => 1, SUFFIX => '.bam', TMPDIR => 1);
    my $bam_pos_sort = Genome::Model::Tools::Sam::SortBam->create(file_name => $self->bam_file, output_file => $sorted_bam);
    $bam_pos_sort->execute() || die "Can't create sorted bam file with Genome::Model::Tools::Sam::SortBam";

    return Bio::DB::Bam->open( $sorted_bam );
}


sub _resolve_output_file {
    my $self = shift;

    $self->debug_message('Resolving output file...');

    unless ( $self->output_directory || $self->output_file) {
        die 'Failed to provide either output_file or output_directory!';
    }

    return $self->output_file if $self->output_file;

    my $resolved_output_file;

    my $output_directory = $self->output_directory;
    unless (-d $output_directory) {
        Genome::Sys->create_directory($output_directory)
            || die "Cannot create output directory '$output_directory'";
    }

    my ($bam_basename,$bam_dirname,$bam_suffix) = File::Basename::fileparse($self->bam_file,qw/\.bam/);
    unless(defined($bam_suffix)) {
        die ('Failed to recognize bam_file '. $self->bam_file .' without bam suffix');
    }
    $resolved_output_file = $output_directory .'/'. $bam_basename;

    my ($bed_basename,$bed_dirname,$bed_suffix) = File::Basename::fileparse($self->bed_file,qw/\.bed/);
    unless(defined($bed_suffix)) {
        die ('Failed to recognize bed_file '. $self->bed_file .' without bed suffix');
    }
    $resolved_output_file .= '-'. $bed_basename;

    if (defined($self->wingspan)) {
        $resolved_output_file .= '-wingspan_'. $self->wingspan;
    }

    $resolved_output_file .= '-alignment_summary-v2.tsv';
    $self->output_file($resolved_output_file);
    return $resolved_output_file;
}

my @output_headers = qw(
    total_bp total_aligned_bp total_unaligned_bp total_duplicate_bp paired_end_bp
    read_1_bp read_2_bp mapped_paired_end_bp proper_paired_end_bp singleton_bp
    total_target_aligned_bp unique_target_aligned_bp duplicate_target_aligned_bp
    total_off_target_aligned_bp unique_off_target_aligned_bp duplicate_off_target_aligned_bp
);

sub _save_summary {
    my ($self, $stats_summary, $output_file) = @_;

    $self->debug_message("Saving summary to output file $output_file");

    my $output_fh = Genome::Sys->open_file_for_writing($self->output_file);
    $output_fh->print(join("\t", @output_headers), "\n");
    no warnings 'uninitialized';
    $output_fh->print(join("\t", @$stats_summary{@output_headers}), "\n");

    return 1;
}


sub _unsorted_bed_file_name_for_chrom {
    my($self,$bed_dir,$chrom) = @_;
    return $bed_dir . "/chr${chrom}.unsorted.bed";
}

sub _sorted_bed_file_name_for_chrom {
    my($self,$bed_dir,$chrom) = @_;
    return $bed_dir . "/chr${chrom}.bed";
}

sub _prepare_split_sorted_bed_files {
    my $self = shift;
    my $wingspan = shift;

    $wingspan ||= 0;

    $self->debug_message('Preparing sorted bed files...');
    # bed file format is
    # chrom start stop name

    # Make a file for each chromosome, and sort by the start and end position 

    my $original = IO::File->new($self->bed_file, 'r') || die "Can't open bed file ".$self->bed_file." for reading: $!";
    my $temp_dir = File::Temp::tempdir( CLEANUP => 1 );

    # First, bin each line by chromosome and make unsorted files
    my %open_handles;
    while(my $line = $original->getline()) {
        $line =~ m/^(\S+)/;  # First column is the chrom name
        my $chrom = $1;

        unless ($open_handles{$chrom}) {
            my $filename = $self->_unsorted_bed_file_name_for_chrom($temp_dir, $chrom);
            $open_handles{$chrom} = IO::File->new($filename, 'w');
            $open_handles{$chrom} || die "Can't open file $filename for writing: $!";
        }

        chomp($line);
        my @line = split(/\s+/, $line);

        # Tack on the wingspan length to each side of the ROI
        $line[1] -= $wingspan;
        $line[2] += $wingspan;
        $open_handles{$chrom}->print(join("\t",@line),"\n");
    }
    $_->close foreach values %open_handles;

    # Now sort each file
    foreach my $chrom ( keys %open_handles ) {
        $self->debug_message("Sorting split bed file for chrom $chrom");

        my $read = IO::File->new($self->_unsorted_bed_file_name_for_chrom($temp_dir, $chrom), 'r');
        $read || die "Can't open file ".$self->_unsorted_bed_file_name_for_chrom($temp_dir, $chrom) . "for reading: $!";

        my @lines;
        while (my $line = $read->getline) {
            chomp $line;
            my @data = split(/\s+/, $line);
            push @lines, \@data;
        }
        $read->close();

        @lines = sort { ($a->[1] <=> $b->[1]) || ($a->[2] <=> $b->[2]) } @lines;
        my $write = IO::File->new($self->_sorted_bed_file_name_for_chrom($temp_dir, $chrom), 'w');
        $write || die "Can't open file ". $self->_sorted_bed_file_name_for_chrom($temp_dir, $chrom) . " for writing: $!";

        foreach ( @lines ) {
            $write->print(join("\t",@$_),"\n");
        }
        $write->close();

        unlink($self->_unsorted_bed_file_name_for_chrom($temp_dir, $chrom));  # Don't need the unsorted file any more
    }

    return $self->{'_sorted_bed_directory'} = $temp_dir;
}


1;
