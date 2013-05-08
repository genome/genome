package Genome::InstrumentData::AlignmentResult::Bwasw;

use strict;
use warnings;
use Carp qw/confess/;
use Data::Dumper;
use File::Basename;
use Genome;
use Getopt::Long;

class Genome::InstrumentData::AlignmentResult::Bwasw {
    is => 'Genome::InstrumentData::AlignmentResult',
    has_constant => [
        aligner_name => { value => 'bwasw', is_param=>1 },
    ],
    has_transient_optional => [
        _bwa_sam_cmd => { is => 'Text' }
    ]
};

sub required_arch_os { 'x86_64' }

sub required_rusage {
    my $class = shift;
    my %p = @_;
    my $instrument_data = delete $p{instrument_data};

    my $tmp_mb = $class->tmp_megabytes_estimated($instrument_data);
    my $mem_mb = 1024 * 16; 
    my $cpus = 4;

    my $mem_kb = $mem_mb*1024;
    my $tmp_gb = $tmp_mb/1024;

    my $user = getpwuid($<);
    my $queue = 'alignment';
    $queue = 'alignment-pd' if (Genome::Config->should_use_alignment_pd);

    my $host_groups;
    $host_groups = qx(bqueues -l $queue | grep ^HOSTS:);
    $host_groups =~ s/\/\s+/\ /;
    $host_groups =~ s/^HOSTS:\s+//;

    my $select  = "select[ncpus >= $cpus && mem >= $mem_mb && gtmp >= $tmp_gb] span[hosts=1]";
    my $rusage  = "rusage[mem=$mem_mb, gtmp=$tmp_gb]";
    my $options = "-M $mem_kb -n $cpus -q $queue";

    my $required_usage = "-R '$select $rusage' $options";

    #check to see if our resource requests are feasible (This uses "maxmem" to check theoretical availability)
    #factor of four is based on current six jobs per host policy this should be revisited later
    my $select_check = "select[ncpus >= $cpus && maxmem >= " . ($mem_mb * 4) . " && maxgtmp >= $tmp_gb] span[hosts=1]";
    my $select_cmd = "bhosts -R '$select_check' $host_groups | grep ^blade";

    my @selected_blades = qx($select_cmd);

    if (@selected_blades) {
        return $required_usage;
    } else {
        die $class->error_message("Failed to find hosts that meet resource requirements ($required_usage). [Looked with `$select_cmd`]");
    }
}

# TODO copied verbatim from normal bwa, but this may be totally off for bwasw
sub tmp_megabytes_estimated {
    my $class = shift || die;
    my $instrument_data = shift;

    my $default_megabytes = 90000;


    if (not defined $instrument_data) {
        return $default_megabytes;
    } elsif ($instrument_data->bam_path) {
        my $bam_path = $instrument_data->bam_path;

        my $scale_factor = 3.25; # assumption: up to 3x during sort/fixmate/sort and also during fastq extraction (2x) + bam = 3

        my $bam_bytes = -s $bam_path;
        unless ($bam_bytes) {
            die $class->error_message("Instrument Data " . $instrument_data->id  . " has BAM ($bam_path) but has no size!");
        }

        if ($instrument_data->can('get_segments')) {
            my $bam_segments = scalar $instrument_data->get_segments;
            if ($bam_segments > 1) {
                $scale_factor = $scale_factor / $bam_segments;
            }
        }

        return int(($bam_bytes * $scale_factor) / 1024**2);
    } elsif ($instrument_data->can("calculate_alignment_estimated_kb_usage")) {
        my $kb_usage = $instrument_data->calculate_alignment_estimated_kb_usage;
        return int(($kb_usage * 3) / 1024) + 100; # assumption: 2x for the quality conversion, one gets rm'di after; aligner generates 0.5 (1.5 total now); rm orig; sort and merge maybe 2-2.5
    } else {
        return $default_megabytes;
    }

    return;
}

# override this from AlignmentResult.pm to filter reads with secondary alignment flag (0x100)
sub _check_read_count {
    my ($self) = @_;
    my $fq_rd_ct = $self->_fastq_read_count;
    my $sam_path = Genome::Model::Tools::Sam->path_for_samtools_version($self->samtools_version);

    my $cmd = "$sam_path view -F 256 -c " . $self->temp_staging_directory . "/all_sequences.bam";
    my $bam_read_count = `$cmd`;
    my $check = "Read count from bam: $bam_read_count and fastq: $fq_rd_ct";

    unless ($fq_rd_ct == $bam_read_count) {
        $self->error_message("$check does not match.");
        return;
    }
    $self->status_message("$check matches.");
    return 1;
}

sub _run_aligner {
    my $self = shift;
    my @input_pathnames = @_;

    # process inputs
    if (@input_pathnames != 1 and @input_pathnames != 2) {
        $self->error_message(
            "Expected 1 or 2 input path names. Got: " . Dumper(\@input_pathnames));
    }
    my $input_filenames = join(' ', @input_pathnames);

    # get temp dir
    my $tmp_dir = $self->temp_scratch_directory;
    my $log_filename = $tmp_dir . '/aligner.log';
    my $raw_sam = $tmp_dir . '/raw_sequences.sam';
    my $fixed_sam = $tmp_dir . '/fixed_sequences.sam';
    my $final_sam = $tmp_dir . '/all_sequences.sam';

    # get refseq info
    my $reference_build = $self->reference_build;
    my $reference_fasta_path = $self->get_reference_sequence_index->full_consensus_path('fa');

    # get params and verify the aligner
    my $params = $self->decomposed_aligner_params;

    my $aligner_version = $self->aligner_version;

    unless (Genome::Model::Tools::Bwa->supports_bwasw($aligner_version)) {
        die $self->error_message("The pipeline does not support using Bwasw with bwa-$aligner_version.");
    }

    my $command_name = Genome::Model::Tools::Bwa->path_for_bwa_version($aligner_version);

    # run cmd
    my $full_command = sprintf '%s bwasw %s %s %s 1>> %s 2>> %s',
        $command_name, $params, $reference_fasta_path, $input_filenames, $raw_sam, $log_filename;

    my $rv = Genome::Sys->shellcmd(
        cmd          => $full_command,
        input_files  => [ @input_pathnames ],
        output_files => [ $raw_sam, $log_filename ],
        skip_if_output_is_present => 0,
    );

    # verify the bwasw logfile
    unless ($self->_verify_bwa_bwasw_did_happen($log_filename)) {
        die $self->error_message("Error running bwasw (unable to verify a successful run of bwasw in the aligner log)");
    }

    if ($rv != 1) {
        die $self->error_message("Error running bwasw (didn't get a good return value from Genome::Sys->shellcmd())");
    }

    my $is_paired = @input_pathnames == 2 ? 1 : 0;
    my $include_secondary = 1;
    my $mark_secondary_as_duplicate = 0;

    $self->status_message("Fixing flags and mates in merged sam file.");

    $self->_fix_sam($raw_sam, $fixed_sam, $is_paired, $include_secondary, $mark_secondary_as_duplicate);

    unlink($raw_sam) || die $self->error_message("Could not unlink $raw_sam.");

    $self->status_message("Resorting fixed sam file by coordinate.");

    my $picard_sort_cmd = Genome::Model::Tools::Picard::SortSam->create(
        sort_order             => 'coordinate',
        input_file             => $fixed_sam,
        output_file            => $final_sam,
        max_records_in_ram     => 2000000,
        maximum_memory         => 8,
        maximum_permgen_memory => 256,
        temp_directory         => $self->temp_scratch_directory,
        use_version            => $self->picard_version,
    );

    # TODO not sure if the following is necessary
    #my $add_rg_cmd = Genome::Model::Tools::Sam::AddReadGroupTag->create(
    #    input_filehandle  => $sorted_sam,
    #    output_filehandle => $final_sam,
    #    read_group_tag    => $self->read_and_platform_group_tag_id,
    #    pass_sam_headers  => 0,
    #);

    unless ($picard_sort_cmd and $picard_sort_cmd->execute) {
        die $self->error_message("Failed to create or execute picard sort command.");
    }

    unlink($fixed_sam) || die $self->error_message("Could not unlink $fixed_sam.");

    return 1;
}


# Main function for fixing sam output.
# TODO This should probably be refactored into smaller functions.
# TODO possible things this overlooks: NM tags, other tags, whether flag 0x2 should be set
# TODO this may be mishandling all the headers in the sam file
sub _fix_sam {
    my $self = shift;
    my $sam_filename = shift;
    my $out_filename = shift;
    my $pe = shift;
    my $include_secondary = shift;
    my $mark_secondary_as_duplicate = shift;

    open my $fh, '<', $sam_filename;
    open my $out_fh, '>', $out_filename;

    my @header_lines = $self->_get_header($fh);
    $self->_write_header($out_fh, @header_lines);

    if ($pe) {
        while (my @read_set = $self->_get_read_set($fh)) {
            die $self->error_message("Too few reads in read set") unless @read_set;
            my $first_primary = shift @read_set;
            my @first_secondary;

            while (@read_set) {
                last if not $self->_equal_seq_qual($first_primary, $read_set[0]);
                push @first_secondary, shift @read_set;
            }

            die $self->error_message("Too few reads in read set") unless @read_set;
            my $last_primary = shift @read_set;
            my @last_secondary;

            while (@read_set) {
                last if not $self->_equal_seq_qual($last_primary, $read_set[0]);
                push @last_secondary, shift @read_set;
            }

            # we should have shifted everything out of @read_set
            die $self->error_message("More than two unique sequences in read set; please make sure bwasw is run without the -H option") if @read_set;

            # verify that the primary reads have the correct flags
            die $self->error_message("Unexpected flags; please make sure bwasw is run with the -M option") if
                ( $self->_get_flag($first_primary, 0x4) and @first_secondary ) or
                ( $self->_get_flag($first_primary, 0x80) ) or
                ( $self->_get_flag($first_primary, 0x100) );

            die $self->error_message("Unexpected flags; please make sure bwasw is run with the -M option") if
                ( $self->_get_flag($last_primary, 0x4) and @last_secondary ) or
                ( $self->_get_flag($last_primary, 0x40) ) or
                ( $self->_get_flag($last_primary, 0x100) );

            # verify that secondary reads have the correct flags
            die $self->error_message("Unexpected flags; please make sure bwasw is run with the -M option") if grep {
                $self->_get_flag($_, 0x4) or
                not $self->_get_flag($_, 0x40) or
                $self->_get_flag($_, 0x80) or
                not $self->_get_flag($_, 0x100)
            } @first_secondary;

            die $self->error_message("Unexpected flags; please make sure bwasw is run with the -M option") if grep {
                $self->_get_flag($_, 0x4) or
                $self->_get_flag($_, 0x40) or
                not $self->_get_flag($_, 0x80) or
                not $self->_get_flag($_, 0x100)
            } @last_secondary;

            # change 0x100 to 0x400 if we want to
            if ($mark_secondary_as_duplicate) {
                for my $read (@first_secondary, @last_secondary) {
                    $self->_unset_flag($read, 0x100);
                    $self->_set_flag($read, 0x400);
                }
            }



            # verify that the primary has an alignment if there are secondary sequences
            die $self->error_message("Primary first sequence was improperly selected by bwasw") if
                @first_secondary and $first_primary->{rname} eq '*';
            die $self->error_message("Primary last sequence was improperly selected by bwasw") if
                @last_secondary and $last_primary->{rname} eq '*';

            # make sure info is properly set
            for my $read ($first_primary, @first_secondary) {
                die $self->error_message("Tlen did not match expected") if
                    $read->{tlen} ne '0' and $read->{tlen} != $self->_calculate_tlen($read, $last_primary, $first_primary);
                # TODO need to set 0x2?
                $self->_set_flag($read, 0x1);
                $self->_set_flag($read, 0x40);
                $self->_add_mate_to_read($read, $last_primary, $first_primary);
            }

            for my $read ($last_primary, @last_secondary) {
                die $self->error_message("Tlen did not match expected") if
                    $read->{tlen} ne '0' and $read->{tlen} != $self->_calculate_tlen($read, $first_primary, $last_primary);
                # TODO need to set 0x2?
                $self->_set_flag($read, 0x1);
                $self->_set_flag($read, 0x80);
                $self->_add_mate_to_read($read, $first_primary, $last_primary);
            }

            if ($include_secondary) {
                $self->_write_read_set($out_fh, ($first_primary, @first_secondary, $last_primary, @last_secondary) );
            } else {
                $self->_write_read_set($out_fh, ($first_primary, $last_primary) );
            }
        }
    } else {
        while (my @read_set = $self->_get_read_set($fh)) {
            die $self->error_message("Too few reads in read set") unless @read_set;
            my $primary = shift @read_set;
            my @secondary;

            while (@read_set) {
                last if not $self->_equal_seq_qual($primary, $read_set[0]);
                push @secondary, shift @read_set;
            }

            # we should have shifted everything out of @read_set
            die $self->error_message("More than two unique sequences in read set; please make sure bwasw is run without the -H option") if @read_set;

            # verify that the primary reads have the correct flags
            die $self->error_message("Unexpected flags; please make sure bwasw is run with the -M option") if
                ( $self->_get_flag($primary, 0x4) and @secondary ) or
                ( $self->_get_flag($primary, 0x100) );

            # verify that secondary reads have the correct flags
            die $self->error_message("Unexpected flags; please make sure bwasw is run with the -M option") if grep {
                $self->_get_flag($_, 0x4) or
                not $self->_get_flag($_, 0x100)
            } @secondary;

            # change 0x100 to 0x400 if we want to
            if ($mark_secondary_as_duplicate) {
                for my $read (@secondary) {
                    $self->_unset_flag($read, 0x100);
                    $self->_set_flag($read, 0x400);
                }
            }

            # verify that the primary has an alignment if there are secondary sequences
            die $self->error_message("Primary sequence was improperly selected by bwasw") if
                @secondary and $primary->{rname} eq '*';

            # make sure that no mate info is set
            for my $read ($primary, @secondary) {
                die $self->error_message("Unexpected flags; flags indicate paired-end data when we were expecting non-paired reads") if
                    $self->_get_flag($read, 0x1) or
                    $self->_get_flag($read, 0x2) or
                    $self->_get_flag($read, 0x8) or
                    $self->_get_flag($read, 0x20) or
                    $self->_get_flag($read, 0x40) or
                    $self->_get_flag($read, 0x80);

                die $self->error_message("Rnext indicates paired-end data when we were expecting non-paired reads") if
                    $read->{rnext} ne '*';

                die $self->error_message("Pnext indicates paired-end data when we were expecting non-paired reads") if
                    $read->{pnext} ne '0';

                die $self->error_message("Tlen indicates paired-end data when we were expecting non-paired reads") if
                    $read->{tlen} ne '0';
            }

            if ($include_secondary) {
                $self->_write_read_set($out_fh, ($primary, @secondary) );
            } else {
                $self->_write_read_set($out_fh, ($primary) );
            }
        }
    }

    close $fh;
    close $out_fh;
}

# set flag to 1 given read and bit
sub _set_flag {
    my $self = shift;
    my $read = shift;
    my $bit = shift;
    my $flag = $read->{flag};
    $flag = $flag | $bit;
    $read->{flag} = $flag;
}

# set flag to 0 given read and bit
sub _unset_flag {
    my $self = shift;
    my $read = shift;
    my $bit = shift;
    my $flag = $read->{flag};
    $flag = $flag & (~$bit);
    $read->{flag} = $flag;
}

# get flag given read and bit
sub _get_flag {
    my $self = shift;
    my $read = shift;
    my $bit = shift;
    return $read->{flag} & $bit;
}

# Note that this takes data from mate and adds it to read, but *not* vice
# versa. Given two reads this must be called twice to fix both records.
sub _add_mate_to_read {
    my $self = shift;
    my $read = shift;
    my $mate = shift;
    my $primary_read = shift;

    if ($mate->{rname} eq $read->{rname} and $read->{rname} ne '*') {
        $read->{rnext} = '=';
    } else {
        $read->{rnext} = $mate->{rname};
    }
    $read->{pnext} = $mate->{pos};

    if ($mate->{flag} & 0x4) {
        $read->{flag} = $read->{flag} | 0x8; # |= wasn't working for some reason
    }

    if ($mate->{flag} & 0x10) {
        $read->{flag} = $read->{flag} | 0x20; # |= wasn't working for some reason
    }
    # TODO verify that bwasw does not reverse complement the read then strip
    # the flag that indicates this

    $read->{tlen} = $self->_calculate_tlen($read, $mate, $primary_read);
}

# calculate tlen given a pair
sub _calculate_tlen {
    my $self = shift;
    my $read = shift;
    my $mate = shift;
    my $primary_read = shift; # this is always the primary alignment, while $read can be primary or secondary

    # Does this seem weird? Yes? Well that's what bwasw does.
    #my $read_len = $self->_cigar_len($read->{cigar});
    my $primary_read_len = $self->_cigar_len($primary_read->{cigar});
    my $mate_len = $self->_cigar_len($mate->{cigar});

    # this is the algorithm bwasw uses for calculating tlen
    if ($read->{rname} eq $mate->{rname}) {
        if ($read->{pnext} + $mate_len > $read->{pos}) {
            return $read->{pnext} + $mate_len - $read->{pos};
        } else {
            return $read->{pnext} - $read->{pos} - $primary_read_len;
        }
    } else {
        return 0;
    }

    # I believe the following block is how TLEN is calculated according to the
    # SAM specification, but this contradicts the behavior of bwasw:
    #if ($read->{rname} eq $mate->{rname}) {
    #    if ($read->{pos} < $mate->{pos}) {
    #        my $length = _cigar_len($read->{cigar}, qw(M D N));
    #        return ($read->{pos} + $length - $mate->{pos} - 1);
    #    } elsif ($read->{pos} > $mate->{pos}) {
    #        my $length = _cigar_len($mate->{cigar}, qw(M D N));
    #        return (-1 * ($mate->{pos} + $length - $read->{pos} - 1));
    #    } else {
    #        # TODO not sure what do in this situation, not sure if this is even valid? which is the "leftmost" read?
    #        return 0;
    #    }
    #} else {
    #    return 0;
    #}
}

# get length of sequence given a cigar string
sub _cigar_len {
    my $self = shift;
    my $cigar = shift;
    my @ops = @_;

    @ops = qw(M D N) unless @ops;

    my $seq_length = 0;

    while ($cigar =~ /([0-9]+)([MIDNSHP=X]+)(.*)/) {
        my $length = $1;
        my $op = $2;
        $cigar = $3;

        $seq_length += $length if grep { $_ eq $op } @ops;
    }

    return $seq_length;
}

# determine if two reads are equal by looking at the sequence and its quality string
sub _equal_seq_qual {
    my $self = shift;
    my $first = shift;
    my $second = shift;

    my $first_seq = $first->{seq};
    my $first_qual = $first->{qual};
    my $second_seq = $second->{seq};
    my $second_qual = $second->{qual};

    if ($self->_get_flag($first, 0x10)) {
        $first_seq = $self->_rev_seq($first->{seq});
        $first_qual = $self->_rev_qual($first->{qual});
    }
    if ($self->_get_flag($second, 0x10)) {
        $second_seq = $self->_rev_seq($second->{seq});
        $second_qual = $self->_rev_qual($second->{qual});
    }

    return $first_seq eq $second_seq and $first_qual eq $second_qual;
}

# reverse a quality string
sub _rev_qual {
    my $self = shift;
    my $seq = shift;
    return join('', reverse(split('', $seq)));
}

# reverse complement a sequence
sub _rev_seq {
    my $self = shift;
    my $seq = shift;

    # TODO is this comprehensive?
    my %complement = (
        'A' => 'T',
        'T' => 'A',
        'G' => 'C',
        'C' => 'G',
    );

    return join('', reverse(map {exists($complement{$_}) ? $complement{$_} : $_ } split('', $seq)));
}

# for sorting by quality
sub _by_quality {
    my $a_unmapped = $a->{flag} & 0x4;
    my $b_unmapped = $b->{flag} & 0x4;
    if ($a_unmapped and $b_unmapped) {
        return 0;
    } elsif ($a_unmapped) {
        return -1;
    } elsif ($b_unmapped) {
        return 1;
    } else {
        return $a->{mapq} <=> $b->{mapq};
    }
}

# Gets the header. Returns an empty list if no line begins with @.
sub _get_header {
    my $self = shift;
    my $fh = shift;

    my @headers;

    my $mark = tell $fh;
    while (my $line = <$fh>) {
        last if $line !~ /^@/;
        chomp $line;
        push @headers, $line;
        $mark = tell $fh;
    }
    seek $fh, $mark, 0;

    return @headers;
}

# Uses _get_read and _get_matching_read to get the next "set" of reads, where a
# set of reads all have the same qname. Returns an empty list if there are no
# reads left, and dies when it encounters a header or a malformed read.
sub _get_read_set {
    my $self = shift;
    my $fh = shift;

    my @read_list;

    my $initial_read = $self->_get_read($fh);

    if ($initial_read) {
        push @read_list, $initial_read;

        while (my $additional_read = $self->_get_matching_read($fh, $initial_read)) {
            push @read_list, $additional_read;
        }
    }

    return @read_list;
}

# Either returns the read as a hashref or undef if there is no read. If it
# encounters a header or malformed read, it will die.
sub _get_read {
    my $self = shift;
    my $fh = shift;

    my $line = <$fh>;

    return undef unless defined $line;

    die $self->error_message("Unexpected header in sam file") if $line =~ /^@/;
    chomp $line;

    my @keys = qw(qname flag rname pos mapq cigar rnext pnext tlen seq qual tags);
    my @fields = split /\t/, $line, 12;

    # we can have an empty tags field, so we expect 11 or 12 fields
    die $self->error_message("Incorrect number of fields in sam line") unless @fields == 11 or @fields == 12;

    my %record;

    $record{$keys[$_]} = $fields[$_] for (0..$#keys);

    $self->_validate_read(\%record);

    return \%record;
}

# Validate a read. This only tests whether a read matches the SAM
# specification. This is not comprehensive. For instance, this does not
# validate flags, optional tags, header, whether things map off the reference,
# etc.
sub _validate_read {
    my $self = shift;
    my $read = shift;

    my $num_regex = '^-?\d+\z';

    # these are just copied out of the sam specification...
    die $self->error_message("Invalid qname") if $read->{qname} !~ /^[!-?A-~]{1,255}\z/;
    die $self->error_message("Invalid flag") if $read->{flag} !~ /$num_regex/;
    die $self->error_message("Invalid flag") if $read->{flag} < 0 or $read->{flag} >= 65536; # 2^16
    die $self->error_message("Invalid rname") if $read->{rname} !~ /(?:\*|[!-()+-<>-~][!-~]*)/;
    die $self->error_message("Invalid pos") if $read->{pos} !~ /$num_regex/;
    die $self->error_message("Invalid pos") if $read->{pos} < 0 or $read->{pos} >= 536870912; # 2^29
    die $self->error_message("Invalid mapq") if $read->{mapq} !~ /$num_regex/;
    die $self->error_message("Invalid mapq") if $read->{mapq} < 0 or $read->{mapq} >= 256; # 2^8
    die $self->error_message("Invalid cigar") if $read->{cigar} !~ /(?:\*|(?:[0-9]+[MIDNSHPX=])+)/;
    die $self->error_message("Invalid rnext") if $read->{rnext} !~ /(?:\*|=|[!-()+-<>-~][!-~]*)/;
    die $self->error_message("Invalid pnext") if $read->{pnext} !~ /$num_regex/;
    die $self->error_message("Invalid pnext") if $read->{pnext} < 0 or $read->{pnext} >= 536870912; # 2^29
    die $self->error_message("Invalid tlen") if $read->{tlen} !~ /$num_regex/;
    die $self->error_message("Invalid tlen") if $read->{tlen} <= -536870912 or $read->{tlen} >= 536870912; # 2^29
    die $self->error_message("Invalid seq") if $read->{seq} !~ /(?:\*|[A-Za-z=.]+)/;
    die $self->error_message("Invalid qual") if $read->{qual} !~ /[!-~]+/;
    # ...better tests here if you want...
}

# Given a read, will return the next read in the file if it has a matching
# qname. Otherwise it will seek back to the original position and return undef.
# Will die if it encounters a header or a malformed read.
sub _get_matching_read {
    my $self = shift;
    my $fh = shift;
    my $cur_read = shift;

    my $mark = tell $fh;

    my $next_read = $self->_get_read($fh);

    if (not defined $next_read) {
        seek $fh, $mark, 0; # reset our position
        return undef;
    }

    my $cur_qname = $cur_read->{qname};
    $cur_qname =~ s/\/[12]$//; # should be unnecessary, but trim it anyways?

    my $next_qname = $next_read->{qname};
    $next_qname =~ s/\/[12]$//; # should be unnecessary, but trim it anyways?

    if ($next_qname ne $cur_qname) {
        seek $fh, $mark, 0; # reset our position
        return undef;
    }

    return $next_read;
}

sub _write_read_set {
    my $self = shift;
    my $fh = shift;
    my @read_set = @_;

    for (@read_set) {
        $self->_write_read($fh, $_);
    }
}

sub _write_read {
    my $self = shift;
    my $fh = shift;
    my $read = shift;

    my @keys = qw(qname flag rname pos mapq cigar rnext pnext tlen seq qual tags);
    my @values;

    for (@keys) {
        if (not defined $read->{$_}) {
            die $self->error_message("Missing field in sam record") if $_ ne 'tags';
        } else {
            push @values, $read->{$_};
        }
    }

    my $sam_line = join "\t", @values;
    print $fh "$sam_line\n";
}

sub _write_header {
    my $self = shift;
    my $fh = shift;
    my @headers = @_;

    for (@headers) {
        print $fh "$_\n";
    }
}

sub _verify_bwa_bwasw_did_happen {
    my ($self, $log_file) = @_;

    unless ($log_file and -e $log_file) {
        $self->error_message("Log file $log_file is does not exist.");
        return;
    }

    unless ($log_file and -s $log_file) {
        $self->error_message("Log file $log_file is empty.");
        return;
    }

    my $line_count = 100;
    my @last_lines = `tail -$line_count $log_file`;

    if (not (
        ($last_lines[-3] =~ /^\[main\] Version:/) and
        ($last_lines[-2] =~ /^\[main\] CMD:/) and
        ($last_lines[-1] =~ /^\[main\] Real time:/) )
    ) {
        $self->error_message("Last lines of $log_file were unexpected. Dumping last $line_count lines.");
        $self->status_message($_) for @last_lines;
        return;
    }
    return 1;
}

sub decomposed_aligner_params {
    my $self = shift;
    my $param_string = $self->aligner_params || '';

    my $param_hash = $self->get_aligner_params_hash($param_string);

    my $cpu_count = $self->_available_cpu_count;
    my $processed_param_string = $self->join_aligner_params_hash($param_hash);

    $self->status_message("[decomposed_aligner_params] cpu count is $cpu_count");
    $self->status_message("[decomposed_aligner_params] bwa bwasw params are: $processed_param_string");

    if ($param_hash->{t} ne $cpu_count) {
        $param_hash->{t} = $cpu_count;
        my $modified_param_string = $self->join_aligner_params_hash($param_hash);
        $self->status_message("[decomposed_aligner_params] autocalculated CPU requirement, bwa bwasw params modified: $modified_param_string");
    }

    if (not exists $param_hash->{M}) {
        $param_hash->{M} = '';
        my $modified_param_string = $self->join_aligner_params_hash($param_hash);
        $self->status_message("[decomposed_aligner_params] forcing -M, bwa bwasw params modified: $modified_param_string");
    }

    my $final_param_string = $self->join_aligner_params_hash($param_hash);

    return $final_param_string;
}

sub aligner_params_for_sam_header {
    my $self = shift;

    my $param_string = $self->aligner_params || '';
    my $param_hash = $self->get_aligner_params_hash($param_string);

    delete $param_hash->{t}; # we don't want cpu count to be part of the sam header

    my $modified_param_string = $self->join_aligner_params_hash($param_hash);

    return "bwa bwasw $modified_param_string";
}

# helper for decomposed_aligner_params and aligner_params_for_sam_header
sub get_aligner_params_hash {
    my $self = shift;
    my $param_string = shift;

    Getopt::Long::Configure("bundling");

    my %param_hash;
    my $rv = Getopt::Long::GetOptionsFromString($param_string,
        'a=i' => \$param_hash{a},
        'b=i' => \$param_hash{b},
        'q=i' => \$param_hash{q},
        'r=i' => \$param_hash{r},
        'w=i' => \$param_hash{w},
        'm=f' => \$param_hash{m},
        't=i' => \$param_hash{t},
        'f=s' => \$param_hash{f},
        'H' => \$param_hash{H},
        'M' => \$param_hash{M},
        'S' => \$param_hash{S},
        'I=i' => \$param_hash{I},
        'T=i' => \$param_hash{T},
        'c=f' => \$param_hash{c},
        'z=i' => \$param_hash{z},
        's=i' => \$param_hash{s},
        'N=i' => \$param_hash{N},
    );

    die $self->error_message("Failed to parse parameter string: $param_string") unless $rv;

    my @switches = qw(H M S);

    for my $key (keys %param_hash) {
        if (not defined $param_hash{$key}) {
            delete $param_hash{$key};
            next;
        }
        if (grep { $key eq $_ } @switches) {
            if ($param_hash{$key} == 1) {
                $param_hash{$key} = '';
            } else {
                delete $param_hash{$key};
            }
        }
    }

    return \%param_hash;
}

# helper for decomposed_aligner_params and aligner_params_for_sam_header
sub join_aligner_params_hash {
    my $self = shift;
    my $param_hash = shift;

    my @param_list;

    for my $key (sort { $a cmp $b } keys %$param_hash) {
        my $val = $param_hash->{$key};
        push @param_list, "-$key";
        push @param_list, $val if $val;
    }

    return join ' ', @param_list;
}

sub fillmd_for_sam {
    return 1;
}

sub requires_read_group_addition {
    return 0;
}

sub supports_streaming_to_bam {
    return 0;
}

sub multiple_reference_mode {
    return 0;
}

sub accepts_bam_input {
    return 0;
}

# Bwasw should just find a corresponding Bwa index and symlink it. This is the
# best we can do within the existing framework if we don't want to recreate an
# identical index already created by the 'regular' bwa module.
sub prepare_reference_sequence_index {
    my $class = shift;
    my $refindex = shift;

    my $staging_dir = $refindex->temp_staging_directory;

    Genome::Sys->create_symlink($refindex->reference_build->get_sequence_dictionary("sam"), $staging_dir ."/all_sequences.dict" );

    my $bwa_index = Genome::Model::Build::ReferenceSequence::AlignerIndex->get_or_create(
        reference_build_id => $refindex->reference_build_id,
        aligner_name       => 'bwa',
        #aligner_params     => $refindex->aligner_params, # none of the aligner params should affect the index step so I think this okay
        aligner_version    => $refindex->aligner_version,
        test_name          => $ENV{GENOME_ALIGNER_INDEX_TEST_NAME},
    );

    for my $filepath (glob($bwa_index->output_dir . "/*")){
        my $filename = File::Basename::fileparse($filepath);
        next if $filename eq 'all_sequences.fa';
        next if $filename eq 'all_sequences.dict';
        Genome::Sys->create_symlink($filepath, $staging_dir . "/$filename");
    }

    $bwa_index->add_user(
        label => 'uses',
        user  => $refindex
    );

    return 1;
}

1;

