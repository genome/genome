package Genome::InstrumentData::AlignmentResult::Bwasw;

use strict;
use warnings;
use Carp qw/confess/;
use Data::Dumper;
use File::Basename;
use File::Copy qw/move/;
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
    my $aligner_params  = delete $p{aligner_params};

    my $tmp_mb = $class->tmp_megabytes_estimated($instrument_data);
    my $mem_mb = 1024 * 16; 
    my $cpus = 4;

    if ($aligner_params and $aligner_params =~ /-t\s*([0-9]+)/) {
        $cpus = $1;
    }

    my $mem_kb = $mem_mb*1024;
    my $tmp_gb = $tmp_mb/1024;

    my $user = getpwuid($<);
    my $queue = $ENV{GENOME_LSF_QUEUE_ALIGNMENT_DEFAULT};
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
    my @input_paths = @_;

    # Process inputs.
    if (@input_paths != 1 and @input_paths != 2) {
        $self->error_message(
            "Expected 1 or 2 input path names. Got: " . Dumper(\@input_paths));
    }
    my $reference_fasta_path = $self->get_reference_sequence_index->full_consensus_path('fa');

    # Get tmp dir.
    my $tmp_dir = $self->temp_scratch_directory;
    my $log_path      = $tmp_dir . '/aligner.log';
    my $raw_sequences = $tmp_dir . '/raw_sequences.sam';
    my $all_sequences = $tmp_dir . '/all_sequences.sam';

    # Verify the aligner and get params.
    my $aligner_version = $self->aligner_version;
    unless (Genome::Model::Tools::Bwa->supports_bwasw($aligner_version)) {
        die $self->error_message(
            "The pipeline does not support using " .
            "bwasw with bwa-$aligner_version."
        );
    }
    my $cmd_path = Genome::Model::Tools::Bwa->path_for_bwa_version($aligner_version);
    my $params = $self->decomposed_aligner_params;

    # Verify inputs and outputs.
    for (@input_paths, $reference_fasta_path) {
        die $self->error_message("Missing input '$_'.") unless -e $_;
        die $self->error_message("Input '$_' is empty.") unless -s $_;
    }

    # Disconnect from db
    $self->_disconnect_from_db();

    # Run bwasw.
    my $full_command = sprintf '%s bwasw %s %s %s 1>> %s 2>> %s',
        $cmd_path, $params, $reference_fasta_path,
        (join ' ', @input_paths), $raw_sequences, $log_path;

    my $rv = Genome::Sys->shellcmd(
        cmd          => $full_command,
        input_files  => [ $reference_fasta_path, @input_paths ],
        output_files => [ $raw_sequences, $log_path ],
        skip_if_output_is_present => 0,
    );

    # Hopefully we're still disconnected
    $self->_check_db_connection();

    # TODO One potential improvement is to stream the output from bwasw
    # straight into _fix_sam. However, _fix_sam currently seeks back and forth
    # while trying to get a read "set" (all primary and secondary alignments
    # for a pair of reads).
    # $self->status_message("Running '$full_command' and streaming output.");
    # my $command_output_fh = IO::File->new( "$full_command |" );

    # Verify the bwasw logfile.
    unless ($self->_verify_bwa_bwasw_did_happen($log_path)) {
        die $self->error_message(
            "Error running bwasw (unable to verify a successful " .
            "run of bwasw in the aligner log)"
        );
    }

    unless ($rv == 1) {
        die $self->error_message(
            "Error running bwasw (didn't get a good return " .
            "value from Genome::Sys->shellcmd)"
        );
    }

    # Fix raw_sequences.sam and write the fixed records to all_sequences.sam.
    $self->status_message("Fixing flags and mates in merged sam file.");

    my $is_paired = @input_paths == 2 ? 1 : 0;
    my $include_secondary = 1;
    my $mark_secondary_as_duplicate = 0;

    $self->_fix_sam($raw_sequences, $all_sequences, $is_paired, $include_secondary, $mark_secondary_as_duplicate);
    unlink($raw_sequences) || $self->status_message("Could not unlink $raw_sequences.");

    # Sort all_sequences.sam.
    $self->status_message("Resorting fixed sam file by coordinate.");
    $self->_sort_sam($all_sequences);

    return 1;
}

sub _disconnect_from_db {
    my ($self) = @_;

    $self->debug_message("Closing data source db handle...");
    if ($self->__meta__->data_source->has_default_handle) {
        if ($self->__meta__->data_source->disconnect_default_handle) {
            $self->debug_message("Disconnected data source db handle (as expected).");
        } else {
            $self->debug_message("Unable to disconnect data source db handle.");
        }
    } else {
        $self->debug_message("Data source db handle already closed.");
    }
}

sub _check_db_connection {
    my ($self) = @_;

    if ($self->__meta__->data_source->has_default_handle) {
        $self->debug_message("Data source db handle unexpectedly reconnected itself.");
    } else {
        $self->debug_message("Data source db handle still closed (as expected).");
    }
}

# Sort a sam file.
sub _sort_sam {
    my ($self, $given_sam) = @_;

    my $unsorted_sam = "$given_sam.unsorted";

    # Prepare sort command
    unless (move($given_sam, $unsorted_sam)) {
        die $self->error_message(
            "Unable to move $given_sam to $unsorted_sam. " .
            "Cannot proceed with sorting.");
    }

    my $picard_sort_cmd = Genome::Model::Tools::Picard::SortSam->create(
        sort_order             => 'coordinate',
        input_file             => $unsorted_sam,
        output_file            => $given_sam,
        max_records_in_ram     => 2000000,
        maximum_memory         => 8,
        maximum_permgen_memory => 256,
        temp_directory         => $self->temp_scratch_directory,
        use_version            => $self->picard_version,
    );

    # Disconnect from db
    $self->_disconnect_from_db();

    # Run sort command
    unless ($picard_sort_cmd and $picard_sort_cmd->execute) {
        die $self->error_message(
            "Failed to create or execute Picard sort command.");
    }

    # Hopefully we're still disconnected
    $self->_check_db_connection();

    # Clean up
    unless (unlink($unsorted_sam)) {
        $self->status_message("Could not unlink $unsorted_sam.");
    }

    return $given_sam;
}

# Main loop for fixing sam output. Opens filehandles, gets headers, pulls in
# reads, modifies them, and then prints them back out to a new filehandle.
sub _fix_sam {
    my ($self, $raw_sequences, $all_sequences, $is_paired, $include_secondary, $mark_secondary_as_duplicate) = @_;

    # Generate RG and PG tag
    my $rg_id = $self->read_and_platform_group_tag_id;
    my $rg_tag = "RG:Z:$rg_id\tPG:Z:$rg_id";

    # Open read filehandles
    open my $all_sequences_read_fh, '<', $all_sequences;
    open my $raw_sequences_read_fh, '<', $raw_sequences;

    # Get existing headers
    my @all_sequences_headers = _get_headers($all_sequences_read_fh); # get headers provided in all_sequences.fa
    my @raw_sequences_headers = _get_headers($raw_sequences_read_fh); # get headers created by Bwa

    # Open write (append) filehandle
    close $all_sequences_read_fh;
    open my $all_sequences_append_fh, '>>', $all_sequences;

    # Create final header
    for my $raw_header (@raw_sequences_headers) {
        unless (grep { index($_, $raw_header) >= 0 } @all_sequences_headers) {
            # Add any unique headers from Bwa to all_sequences.fa
            $self->status_message("Didn't find '$raw_header' in all_sequences.sam. Adding it. This probably isn't a problem.");
            print $all_sequences_append_fh "$raw_header\n";
        }
    }

    # Disconnect from db before processing the rest of the sam file
    $self->_disconnect_from_db();

    # Loop over all records
    while (my @read_set = _get_read_set($raw_sequences_read_fh)) {
        die "No reads in read set; this should not happen"
            unless @read_set;
        my $read_set_dump = Dumper(\@read_set);

        eval {
            my @fixed_read_set = $is_paired
                ? _fix_paired_set($include_secondary, $mark_secondary_as_duplicate, @read_set)
                : _fix_unpaired_set($include_secondary, $mark_secondary_as_duplicate, @read_set);

            _add_rg_tag($_, $rg_tag) for @fixed_read_set;
            print $all_sequences_append_fh _read_to_string($_) for @fixed_read_set;
        };
        if ($@) {
            my $error = sprintf "Error fixing SAM: %s\nCurrent read set: %s",
                $@, $read_set_dump;
            die $self->error_message($error);
        }
    }

    close $raw_sequences_read_fh;
    close $all_sequences_append_fh;

    # Hopefully we're still disconnected
    $self->_check_db_connection();
}

# Fixes a read set containing paired reads.
# TODO possible things this overlooks: NM tags, whether flag 0x2 should be set
sub _fix_paired_set {
    my ($include_secondary, $mark_secondary_as_duplicate, @read_set) = @_;

    # split up the records in our read set
    my $first_primary = shift @read_set;
    my @first_secondary;

    while (@read_set) {
        last if _get_flag($read_set[0], 0x4) or not _get_flag($read_set[0], 0x40);
        push @first_secondary, shift @read_set;
    }

    # mates must still be in @read_set
    die "Too few reads in read set (alignment ended early?)."
        unless @read_set;

    my $last_primary = shift @read_set;
    my @last_secondary;

    while (@read_set) {
        last if _get_flag($read_set[0], 0x4) or not _get_flag($read_set[0], 0x80);
        push @last_secondary, shift @read_set;
    }

    # we should have shifted everything out of @read_set
    die "Too many reads in read set (make sure bwasw is not used with -H)."
        if @read_set;

    # verify that primary alignments have expected flags
    _verify_primary($first_primary, scalar(@first_secondary), [0x1, 0x40], [0x80, 0x100]);
    _verify_primary($last_primary, scalar(@last_secondary), [0x1, 0x80], [0x40, 0x100]);

    # verify that secondary alignments have expected flags, and change
    # 0x100 to 0x400 if we want to on first secondary alignments
    for my $read (@first_secondary) {
        _verify_secondary($read, [0x1, 0x40, 0x100], [0x4, 0x80]);
        _mark_as_duplicate($read) if $mark_secondary_as_duplicate;
    }

    for my $read (@last_secondary) {
        _verify_secondary($read, [0x1, 0x80, 0x100], [0x4, 0x40]);
        _mark_as_duplicate($read) if $mark_secondary_as_duplicate;
    }

    # fill in mate info
    for my $read ($first_primary, @first_secondary) {
        _flag_unmapped($read, [0x1, 0x40]); # TODO need to set 0x2?
        _resolve_mate_info($read, $last_primary, $first_primary);
    }

    for my $read ($last_primary, @last_secondary) {
        _flag_unmapped($read, [0x1, 0x80]); # TODO need to set 0x2?
        _resolve_mate_info($read, $first_primary, $last_primary);
    }

    # write revised records
    return $include_secondary
        ? ($first_primary, @first_secondary, $last_primary, @last_secondary)
        : ($first_primary, $last_primary);
}

# Fixes a read set containing non-paired reads.
# TODO possible things this overlooks: NM tags, whether flag 0x2 should be set
sub _fix_unpaired_set {
    my ($include_secondary, $mark_secondary_as_duplicate, @read_set) = @_;

    # split up the records in our read set
    my $primary = shift @read_set;
    my @secondary = @read_set;

    # verify that primary alignment has expected flags and no mate info
    _verify_primary($primary, scalar(@secondary), [], [0x1, 0x2, 0x8, 0x20, 0x40, 0x80, 0x100]);
    _verify_no_mate($primary);

    for my $read (@secondary) {
        # verify that secondary alignments have expected flags and no mate info
        _verify_secondary($read, [0x100], [0x1, 0x2, 0x4, 0x8, 0x20, 0x40, 0x80]);
        _verify_no_mate($read);

        # change 0x100 to 0x400 if we want to on secondary alignments
        _mark_as_duplicate($read) if $mark_secondary_as_duplicate;
    }

    # write revised records
    return $include_secondary ? ($primary, @secondary) : $primary;
}

# Note that this takes data from mate and adds it to read, but *not* vice
# versa. Given two reads this must be called twice to fix both records.
sub _add_mate_to_read {
    my ($read, $mate_primary_read) = @_;

    if ($mate_primary_read->{rname} eq $read->{rname} and $read->{rname} ne '*') {
        $read->{rnext} = '=';
    } else {
        $read->{rnext} = $mate_primary_read->{rname};
    }
    $read->{pnext} = $mate_primary_read->{pos};

    if ($mate_primary_read->{flag} & 0x4) {
        $read->{flag} = $read->{flag} | 0x8; # |= wasn't working for some reason
    }

    if ($mate_primary_read->{flag} & 0x10) {
        $read->{flag} = $read->{flag} | 0x20; # |= wasn't working for some reason
    }
    # TODO verify that bwasw does not reverse complement the read then strip
    # the flag that indicates this
}

# Set the given bit to 1 in a read's flag.
sub _set_flag {
    my ($read, $bit) = @_;
    my $flag = $read->{flag};
    $flag = $flag | $bit;
    $read->{flag} = $flag;
}

# Set the given bit to 0 in a read's flag.
sub _unset_flag {
    my ($read, $bit) = @_;
    my $flag = $read->{flag};
    $flag = $flag & (~$bit);
    $read->{flag} = $flag;
}

# Get whether a given bit is set in a read's flag.
sub _get_flag {
    my ($read, $bit) = @_;
    return $read->{flag} & $bit;
}

# Helper for printing descriptive errors messages about incorrect flags.
sub _flag_description {
    my ($flag) = @_;
    my %description = (
        0x1   => 'template having multiple segments in sequencing',
        0x2   => 'each segment properly aligned according to the aligner',
        0x4   => 'segment unmapped',
        0x8   => 'next segment in the template unmapped',
        0x10  => 'SEQ being reverse complemented',
        0x20  => 'SEQ of the next segment in the template being reversed',
        0x40  => 'the first segment in the template',
        0x80  => 'the last segment in the template',
        0x100 => 'secondary alignment',
        0x200 => 'not passing quality controls',
        0x400 => 'PCR or optical duplicate',
    );
    return $description{$flag};
}

# Verify a primary alignment
sub _verify_primary {
    my ($read, $secondary_count, $good_flag_ref, $bad_flag_ref) = @_;

    my @errors;
    if (_get_flag($read, 0x4)) {
        push @errors, "Primary alignment unmapped but contained flags besides 0x4"
            if $read->{flag} ne '4';
        push @errors, "Primary alignment unmapped but still found secondary alignments"
            if $secondary_count > 0;
    } else {
        push @errors,
            map { sprintf "Primary alignment is missing flag %s (%s).", $_, _flag_description($_) }
            grep { not _get_flag($read, $_) } @$good_flag_ref;
        push @errors,
            map { sprintf "Primary alignment has incorrect flag %s (%s).", $_, _flag_description($_) }
            grep { _get_flag($read, $_) } @$bad_flag_ref;
        push @errors, "Primary alignment is mapped but does not have a proper rname"
            if $read->{rname} eq '*';
    }
    die join ' ', @errors if @errors;
}

# Verify a secondary alignment
sub _verify_secondary {
    my ($read, $good_flag_ref, $bad_flag_ref) = @_;

    my @errors;
    push @errors,
        map { sprintf "Secondary alignment is missing flag %s (%s).", $_, _flag_description($_) }
        grep { not _get_flag($read, $_) } @$good_flag_ref;
    push @errors,
        map { sprintf "Secondary alignment has incorrect flag %s (%s).", $_, _flag_description($_) }
        grep { _get_flag($read, $_) } @$bad_flag_ref;

    die join ' ', @errors if @errors;
}

# Verify a non-paired read does not contain mate information
sub _verify_no_mate {
    my ($read, $read_set_dump) = @_;

    my @errors;
    push @errors, "Rnext indicates paired-end data when we were expecting non-paired reads" if
        $read->{rnext} ne '*';
    push @errors, "Pnext indicates paired-end data when we were expecting non-paired reads" if
        $read->{pnext} ne '0';
    push @errors, "Tlen indicates paired-end data when we were expecting non-paired reads" if
        $read->{tlen} ne '0';

    die join ' ', @errors if @errors;
}

# Change the secondary alignment flag to PCR or optical duplicate flag.
sub _mark_as_duplicate {
    my ($read) = @_;

    _unset_flag($read, 0x100);
    _set_flag($read, 0x400);
}

# Add additional flags to an unmapped read.
sub _flag_unmapped {
    my ($read, $flag_ref) = @_;

    # already verified that reads with 0x4 set have only 0x4 set in _verify_primary
    if (_get_flag($read, 0x4)) {
        for my $flag (@$flag_ref) {
            _set_flag($read, $flag)
        }
    }
}

# Make sure the tlen for paired reads is correct. Fill in the mate info for
# secondary reads.
sub _resolve_mate_info {
    my ($read, $mate_primary_read, $primary_read) = @_;

    my $calculated_tlen = _calculate_tlen($read, $mate_primary_read, $primary_read);

    if ($read->{tlen} ne '0' and $read->{tlen} != $calculated_tlen) {
        die sprintf "Read's tlen did not match calculated tlen: %s != %s",
            $read->{tlen}, $calculated_tlen;
    }

    # This takes a read and fills in the mate information based on the primary
    # mate read.
    _add_mate_to_read($read, $mate_primary_read);

    # Recalculate the tlen, it might have changed after filling in the mate
    # information. (TODO refactor so we're not calculating the tlen twice.)
    $read->{tlen} = _calculate_tlen($read, $mate_primary_read, $primary_read);
}

# calculate tlen given a pair
sub _calculate_tlen {
    my ($read, $mate_primary_read, $primary_read) = @_; # this is always the primary alignment, while $read can be primary or secondary

    # Does this seem weird? Yes? Well that's what bwasw does.
    #my $read_len = _cigar_len($read->{cigar});
    my $primary_read_len = _cigar_len($primary_read->{cigar});
    my $mate_len = _cigar_len($mate_primary_read->{cigar});

    # this is the algorithm bwasw uses for calculating tlen
    if ($read->{rname} eq $mate_primary_read->{rname}) {
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
    #if ($read->{rname} eq $mate_primary_read->{rname}) {
    #    if ($read->{pos} < $mate_primary_read->{pos}) {
    #        my $length = _cigar_len($read->{cigar}, qw(M D N));
    #        return ($read->{pos} + $length - $mate_primary_read->{pos} - 1);
    #    } elsif ($read->{pos} > $mate_primary_read->{pos}) {
    #        my $length = _cigar_len($mate_primary_read->{cigar}, qw(M D N));
    #        return (-1 * ($mate_primary_read->{pos} + $length - $read->{pos} - 1));
    #    } else {
    #        # TODO not sure what do in this situation, not sure if this is even valid? which is the "leftmost" read?
    #        return 0;
    #    }
    #} else {
    #    return 0;
    #}
}

# Determine the length of a sequence for a given cigar string.
sub _cigar_len {
    my ($cigar, @ops) = @_;

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

# Gets the header. Returns an empty list if no line begins with @.
sub _get_headers {
    my ($fh) = @_;
    my @headers;

    my $mark = tell $fh;
    while (my $line = <$fh>) {
        last if $line !~ /^@/;
        push @headers, $line;
        $mark = tell $fh;
    }
    seek $fh, $mark, 0; # rewind to $mark
    chomp for @headers;

    return @headers;
}

# Uses _get_read and _get_matching_read to get the next "set" of reads, where a
# set of reads all have the same qname. Returns an empty list if there are no
# reads left, and dies when it encounters a header or a malformed read.
sub _get_read_set {
    my ($fh) = @_;

    my @read_list;
    my $initial_read = _get_read($fh);

    if ($initial_read) {
        push @read_list, $initial_read;

        while (my $additional_read = _get_matching_read($fh, $initial_read)) {
            push @read_list, $additional_read;
        }
    }

    return @read_list;
}

# Given a read, will return the next read in the file if it has a matching
# qname. Otherwise it will seek back to the original position and return undef.
# Will die if it encounters a header or a malformed read.
sub _get_matching_read {
    my ($fh, $cur_read) = @_;

    my $mark = tell $fh;

    my $next_read = _get_read($fh);

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

# Either returns the read as a hashref or undef if there is no read. If it
# encounters a header or malformed read, it will die.
sub _get_read {
    my ($fh) = @_;

    my $line = <$fh>;

    return undef unless defined $line;

    die "Unexpected header in sam file" if $line =~ /^@/;
    chomp $line;

    my @keys = qw(qname flag rname pos mapq cigar rnext pnext tlen seq qual tags);
    my @fields = split /\t/, $line, 12;

    # we can have an empty tags field, so we expect 11 or 12 fields
    die "Incorrect number of fields in sam line" unless @fields == 11 or @fields == 12;

    my %record;

    $record{$keys[$_]} = $fields[$_] for (0..$#keys);

    _validate_read(\%record);

    return \%record;
}

# Validate a read. This only tests whether a read matches the SAM
# specification. This is not comprehensive. For instance, this does not
# validate flags, optional tags, header, whether things map off the reference,
# etc.
sub _validate_read {
    my ($read) = @_;

    my $num_regex = '^-?\d+\z';

    my @errors;

    # these are just copied out of the sam specification...
    push @errors, "Invalid qname" if $read->{qname} !~ /^[!-?A-~]{1,255}\z/;
    push @errors, "Invalid flag"  if $read->{flag}  !~ /$num_regex/;
    push @errors, "Invalid flag"  if $read->{flag}  <  0
                                  or $read->{flag}  >= 65536; # 2^16
    push @errors, "Invalid rname" if $read->{rname} !~ /(?:\*|[!-()+-<>-~][!-~]*)/;
    push @errors, "Invalid pos"   if $read->{pos}   !~ /$num_regex/;
    push @errors, "Invalid pos"   if $read->{pos}   <  0
                                  or $read->{pos}   >= 536870912; # 2^29
    push @errors, "Invalid mapq"  if $read->{mapq}  !~ /$num_regex/;
    push @errors, "Invalid mapq"  if $read->{mapq}  <  0
                                  or $read->{mapq}  >= 256; # 2^8
    push @errors, "Invalid cigar" if $read->{cigar} !~ /(?:\*|(?:[0-9]+[MIDNSHPX=])+)/;
    push @errors, "Invalid rnext" if $read->{rnext} !~ /(?:\*|=|[!-()+-<>-~][!-~]*)/;
    push @errors, "Invalid pnext" if $read->{pnext} !~ /$num_regex/;
    push @errors, "Invalid pnext" if $read->{pnext} <  0
                                  or $read->{pnext} >= 536870912; # 2^29
    push @errors, "Invalid tlen"  if $read->{tlen}  !~ /$num_regex/;
    push @errors, "Invalid tlen"  if $read->{tlen}  <= -536870912
                                  or $read->{tlen}  >= 536870912; # 2^29
    push @errors, "Invalid seq"   if $read->{seq}   !~ /(?:\*|[A-Za-z=.]+)/;
    push @errors, "Invalid qual"  if $read->{qual}  !~ /[!-~]+/;
    # ...better tests here if you want...

    die join ' ', @errors if @errors;
}

# Add read and project group to tags.
sub _add_rg_tag {
    my ($read, $rg_tag) = @_;
    my $current_tags = $read->{tags};
    $read->{tags} = $current_tags ? "$rg_tag\t$current_tags" : $rg_tag;
}

# Turn a read hash back into a record we can print to a filehandle.
sub _read_to_string {
    my ($read) = @_;

    my @keys = qw(qname flag rname pos mapq cigar rnext pnext tlen seq qual tags);
    my @values;

    for (@keys) {
        if (not defined $read->{$_}) {
            die "Missing field in sam record" if $_ ne 'tags';
        } else {
            push @values, $read->{$_};
        }
    }

    return (join "\t", @values) . "\n";
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

    # Make sure the thread count argument matches the number of CPUs available.
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

