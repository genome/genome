package Genome::Model::Tools::Sam::R2;

use strict;
use warnings;

use Genome;
use Command;
use IO::File;
use POSIX;

class Genome::Model::Tools::Sam::R2{
    is  => 'Genome::Model::Tools::Sam',
    has => [
        input_files => {
            is => 'Text',
            is_many => 1,
            doc => 'SAM or BAM input files to merge.  File extension controls format - a .bam is expected to be a BAM file.  Files with no extension are assumed to be SAM.'
        },
        output_file => {
            is => 'Text',
            doc => 'Merged SAM or BAM file to write.  If the filename ends with .bam, it is written in BAM format.  SAM format is used otherwise.  Use - to write to STDOUT.'
        },
        rand_seed => {
            is_optional => 1,
            is => 'Number',
            doc => 'Random number generator seed.  When multiple alignments or alignment pairs for the same read exist and have the same quality, one is randomly selected.  ' .
                   'If this random selection must be reproducable, specify a value for this argument.'
        }
    ],
};

sub help_brief {
    return 'Tool to merge BAM or SAM files aligned against a split reference (eg one that was larger than 4GiB before being broken up).';
}

sub help_detail {
    return 'Tool to merge BAM or SAM files aligned against a split reference (eg one that was larger than 4GiB before being broken up).  ' .
           'The alignments in the input files must be sorted by read name (samtools sort -n will do this).';
}

sub get_int {
    my $str = shift;
    my $str_idx = shift;
    my $str_len = length($$str);
    my $str_idx_start = $$str_idx;
    for(; $$str_idx < $str_len && isdigit(substr($$str, $$str_idx, 1)); ++$$str_idx) {}
    if($$str_idx == $str_idx_start) {
        die;
    }
    return int(substr($$str, $str_idx_start, $$str_idx - $str_idx_start));
}

sub strnum_cmp
{
    my $lhs = shift;
    my $rhs = shift;

    my $lhs_idx = 0;
    my $rhs_idx = 0;
    my $lhs_end = length($lhs) - 1;
    my $rhs_end = length($rhs) - 1;
    my $lhs_num;
    my $rhs_num;
    my $rhs_curr;
    my $lhs_curr;

    while($lhs_idx < $lhs_end && $rhs_idx < $rhs_end) {
        $lhs_curr = substr($lhs, $lhs_idx, 1);
        $rhs_curr = substr($rhs, $rhs_idx, 1);
        if(isdigit($lhs_curr) && isdigit($rhs_curr)) {
            $lhs_num = get_int(\$lhs, \$lhs_idx);
            $rhs_num = get_int(\$rhs, \$rhs_idx);
            if($lhs_num < $rhs_num) {
                return -1;
            }
            elsif($lhs_num > $rhs_num) {
                return 1;
            }
        }
        else {
            if($lhs_curr ne $rhs_curr) {
                return $lhs_curr cmp $rhs_curr;
            }
            ++$lhs_idx;
            ++$rhs_idx;
        }
    }
    return $lhs_idx <=> $rhs_idx;
}

sub alignment_cmp {
    my $lhs = shift;
    my $rhs = shift;
    my $lhs_mapped = ($lhs->{flag} & 4) == 0;
    my $rhs_mapped = ($rhs->{flag} & 4) == 0;
    if(!$lhs_mapped && $rhs_mapped) {
        return -1;
    }
    if($lhs_mapped && !$rhs_mapped) {
        return 1;
    }
    if(!$lhs_mapped && !$rhs_mapped) {
        return 0;
    }
    return $rhs->{NM} <=> $lhs->{NM};
}

sub alignment_pair_cmp {
    my $lhs = shift;
    my $rhs = shift;
    my $lhs_num_aln = 0;
    my $rhs_num_aln = 0;
    my $lhs_nm = 0;
    my $rhs_nm = 0;
    unless($lhs->[0]{flag} & 4) {
        ++$lhs_num_aln;
        $lhs_nm += $lhs->[0]{NM};
    }
    unless($lhs->[1]{flag} & 4) {
        ++$lhs_num_aln;
        $lhs_nm += $lhs->[1]{NM};
    }
    unless($rhs->[0]{flag} & 4) {
        ++$rhs_num_aln;
        $rhs_nm += $rhs->[0]{NM};
    }
    unless($rhs->[1]{flag} & 4) {
        ++$rhs_num_aln;
        $rhs_nm += $rhs->[1]{NM};
    }
    if($lhs_num_aln != $rhs_num_aln) {
        # One alignment pair has more reads mapped than the other; prefer the one with more
        return $lhs_num_aln <=> $rhs_num_aln;
    }
    if($lhs_num_aln == 2) {
        # Both reads of both alignment pairs aligned
        my $lhs_same_chrom = $lhs->[0]{ref_name} eq $lhs->[1]{ref_name};
        my $rhs_same_chrom = $rhs->[0]{ref_name} eq $rhs->[1]{ref_name};
        if(!$lhs_same_chrom && $rhs_same_chrom) {
            # The reads of lhs mapped to different chromosomes but those of rhs mapped to the same; lhs is considered inferior
            return 1;
        }
        if($lhs_same_chrom && !$rhs_same_chrom) {
            # The reads of lhs mapped to the same chromosome but those of rhs mapped to different ones; rhs is considered inferior
            return -1;
        }
        # For each pair, boths reads mapped to the same chromosome; lhs is considered inferior if it has more total mismatches than rhs
        return $rhs->[0]{NM} + $rhs->[1]{NM} <=> $lhs->[0]{NM} + $lhs->[1]{NM};
    }
    if($lhs_num_aln == 0) {
        # Neither alignment pair has a mapped read; the alignment pairs are equally bad
        return 0;
    }
}

sub open_bamsam_in {
    my $file_struct = shift;
    if($file_struct->{type} eq 'BAM') {
        $file_struct->{file} = new IO::File;
        $file_struct->{file}->open('samtools view -h "' . $file_struct->{name} . '" |');
    }
    elsif($file_struct->{type} eq 'SAM') {
        $file_struct->{file} = new IO::File($file_struct->{name});
    }
    else {
        die 'Unknown type specified for "' . $file_struct->{name} . "\".\n";
    }
    unless($file_struct->{file}) {
        die 'Failed to open "' . $file_struct->{name} . "\"\n.";
    }
    $file_struct->{curr_alignment_index} = -1;
    $file_struct->{curr_alignment} = undef;
    $file_struct->{prev_alignment} = undef;
    $file_struct->{at_eof} = 0;
}

sub open_bamsam_out {
    my $file_struct = shift;
    if($file_struct->{type} eq 'BAM') {
        $file_struct->{file} = new IO::File;
        $file_struct->{file}->open('| samtools view -S -b /dev/stdin > "' . $file_struct->{name} . '"');
    }
    elsif($file_struct->{type} eq 'SAM') {
        $file_struct->{file} = new IO::File($file_struct->{name} eq '-' ? stdout : '> ' . $file_struct->{name});
    }
    else {
        die 'Unknown type specified for "' . $file_struct->{name} . "\".\n";
    }
    unless($file_struct->{file}) {
        die 'Failed to open "' . $file_struct->{name} . "\"\n.";
    }
}

sub get_alignment {
    my $file_struct = shift;
    my $line;
    # Read one alignment or set eof indicator, looping past header and whitespace lines
    for(;;) {
        if(exists $file_struct->{buffline}) {
            $line = $file_struct->{buffline};
            delete $file_struct->{buffline};
        }
        else {
            $line = readline($file_struct->{file});
        }
        if(!defined($line)) {
            $file_struct->{at_eof} = 1;
            last;
        }
        elsif($line !~ /^@/ && $line !~ /^\s*$/) {
            ++$file_struct->{curr_alignment_index};
            chomp $line;
            my @s = split /\t/, $line;
            if(scalar(@s) < 10) {
                die "Failed to parse alignment record " . $file_struct->{curr_alignment_index} . " from \"" . $file_struct->{name} . "\".\n";
            }
            $file_struct->{prev_alignment} = $file_struct->{curr_alignment};
            $file_struct->{curr_alignment} = {};
            $file_struct->{curr_alignment}{file} = $file_struct;
            $file_struct->{curr_alignment}{raw} = $line;
            $file_struct->{curr_alignment}{name} = $s[0];
            $file_struct->{curr_alignment}{flag} = $s[1];
            $file_struct->{curr_alignment}{ref_name} = $s[2];
            $file_struct->{curr_alignment}{ref_pos} = $s[3];
            foreach my $tag (@s[4..$#s]) {
                if($tag =~ /^\s*NM:i:(\d+)\s*$/) {
                    $file_struct->{curr_alignment}{NM} = $1;
                }
            }
            # If the alignment is mapped, verify that it has an NM tag
            unless($file_struct->{curr_alignment}{flag} & 4) {
                unless(exists($file_struct->{curr_alignment}{NM})) {
                    die 'Alignment index ' . $file_struct->{curr_alignment_index} . ' in "' . $file_struct->{name} . "\" is mapped but lacks an NM tag.";
                }
            }
            # Verify that the alignment is in the correct order with respect to the previous alignment
            if( defined($file_struct->{prev_alignment}) && strnum_cmp($file_struct->{curr_alignment}{name}, $file_struct->{prev_alignment}{name}) == -1 ) {
                die 'Alignment index ' . $file_struct->{curr_alignment_index} . ' in "' . $file_struct->{name} . "\" is not properly ordered by read name.\n";
            }
            last;
        }
    }
}

sub execute {
    my $self = shift;

    $self->dump_status_messages(1);
    $self->dump_error_messages(1);

    my @in_filenames = $self->input_files;
    my $out_filename = $self->output_file;

    if(scalar(@in_filenames) < 2) {
        die "Please supply at least two input files.\n";
    }

    if(defined($self->rand_seed)) {
        srand($self->rand_seed);
    }

    # Open input files and output file

    my @in_files;
    foreach my $in_filename (@in_filenames) {
        my $in_file = {file => undef, name => $in_filename, type => ($in_filename =~ /\.bam\s*$/i ? 'BAM' : 'SAM'), at_eof => 0};
        open_bamsam_in($in_file);
        push @in_files, $in_file;
    }

    my %out_file = (file => undef, name => $out_filename, type => ($out_filename =~ /\.bam\s*$/i ? 'BAM' : 'SAM'));
    open_bamsam_out(\%out_file);

    # Merge all headers from every input file

    foreach my $in_file (@in_files) {
        for(;;) {
            my $line = $in_file->{file}->getline();
            if(!defined($line)) {
                $in_file->{at_eof} = 1;
                last;
            }
            if(substr($line, 0, 1) eq '@') {
                $out_file{file}->print($line);
            }
            else {
                $in_file->{buffline} = $line;
                last;
            }
        }
    }

    # Merge alignments

    my @a_heap;
    my @a_p_pool;
    my @a_f_pool;
    my $curr;
    my $partial;
    my $file;
    my $same_end;
    my $same_count;
    my $pool_end;
    my $selected;

    foreach my $in_file (@in_files) {
        get_alignment($in_file);
        unless($in_file->{at_eof}) {
            push @a_heap, $in_file->{curr_alignment};
        }
    }
    @a_heap = sort { strnum_cmp($a->{name}, $b->{name}) } @a_heap;
    for(;;) {
        if(scalar(@a_heap) > 0) {
            $curr = $a_heap[0];
            if(scalar(@a_f_pool) == 0 && scalar(@a_p_pool) == 0 ||
               scalar(@a_f_pool) != 0 && $curr->{name} eq $a_f_pool[0]{name} ||
               scalar(@a_p_pool) != 0 && $curr->{name} eq $a_p_pool[0][0]{name})
            {
                shift @a_heap;
                $file = $curr->{file};
                for(;;) {
                    if(($curr->{flag} & 1) == 0) {
                        push @a_f_pool, $curr;
                    }
                    else {
                        if(defined($partial)) {
                            if( ( ($curr->{flag} ^ $partial->{flag}) & 192 ) != 192 ) {
                                die "Read unmatched pair from from \"" . $file->{name} . "\": pair does not consist of a read1 and read2 (eg it may be two read1s or two read2s).\n";
                            }
                            push @a_p_pool, [$partial, $curr];
                            $partial = undef;
                        }
                        else {
                            # Found half of an alignment pair.  The other half will be discovered in another iteration of this inner
                            # loop (or an error will be thrown).
                            $partial = $curr;
                        }
                    }
                    get_alignment($file);
                    if($file->{at_eof}) {
                        if(defined($partial)) {
                            die "Failed to find second half of paired-end alignment in \"" . $file->{name} . "\": EOF reached.";
                        }
                        last;
                    }
                    if($file->{curr_alignment}{name} ne $curr->{name}) {
                        if(defined($partial)) {
                            die "Failed to find second half of paired-end alignment in \"" . $file->{name} . "\": encountered alignment " .
                                "of a read with a different name.\n";
                        }
                        push @a_heap, $file->{curr_alignment};
                        @a_heap = sort { strnum_cmp($a->{name}, $b->{name}) } @a_heap;
                        last;
                    }
                    $curr = $file->{curr_alignment};
                }
                next;
            }
        }
        if(scalar(@a_p_pool) == 0 && scalar(@a_f_pool) == 0) {
            last;
        }
        unless(scalar(@a_p_pool) == 0) {
            # Sort pool in descending order of pair quality
            @a_p_pool = sort {alignment_pair_cmp($b, $a)} @a_p_pool;
            # Find the first alignment pair in the pool with a different pair quality.  This is the pair one after the
            # last pair with the same mapping quality - we are making the interval [sameStart, sameEnd) containing
            # the alignments with identical best pair quality for the same read pair.
            $same_end = 1;
            $pool_end = scalar(@a_p_pool);
            for(; $same_end < $pool_end && alignment_pair_cmp($a_p_pool[$same_end], $a_p_pool[0]) == 0; ++$same_end) {}
            $same_count = $same_end;
            if($same_count == 1) {
                # One alignment pair has a higher pair quality than any others for the current read pair
                $selected = 0;
            }
            else {
                # Randomly select from the set of alignment pairs tied for best for the current read pair
                $selected = int(rand($same_count + 1));
            }
            $out_file{file}->print($a_p_pool[$selected][0]{raw}, "\n");
            $out_file{file}->print($a_p_pool[$selected][1]{raw}, "\n");
            @a_p_pool = ();
        }
        unless(scalar(@a_f_pool) == 0) {
            # Sort pool in descending order of pair quality
            @a_f_pool = sort {alignment_cmp($b, $a)} @a_f_pool;
            # Find the first alignment in the pool with a different pair quality.  This is the pair one after the
            # last pair with the same mapping quality - we are making the interval [sameStart, sameEnd) containing
            # the alignments with identical best pair quality for the same read pair.
            $same_end = 1;
            $pool_end = scalar(@a_f_pool);
            for(; $same_end < $pool_end && alignment_cmp($a_f_pool[$same_end], $a_f_pool[0]) == 0; ++$same_end) {}
            $same_count = $same_end;
            if($same_count == 1) {
                # One alignment has a higher pair quality than any others for the current read
                $selected = 0;
            }
            else {
                # Randomly select from the set of alignments tied for best for the current read
                $selected = int(rand($same_count + 1));
            }
            $out_file{file}->print($a_f_pool[$selected]{raw}, "\n");
            @a_f_pool = ();
        }
    }
}

1;
