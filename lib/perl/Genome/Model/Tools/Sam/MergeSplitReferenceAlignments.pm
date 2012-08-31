package Genome::Model::Tools::Sam::MergeSplitReferenceAlignments;

use strict;
use warnings;

use Genome;
use Command;
use IO::File;

class Genome::Model::Tools::Sam::MergeSplitReferenceAlignments {
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
        softclip_mismatch => {
            is => 'Int',
            is_optional => 1,
            doc => 'If this flag is specified, softclipped bases will be added to the mismatch count from the NM: tag when choosing the correct read pair',
        },
        rand_seed => {
            is_optional => 1,
            default_value => 12345,
            is => 'Number',
            doc => 'Random number generator seed.  When multiple alignments or alignment pairs for the same read exist and have the same quality, one is randomly selected.  ' .
                   'If this random selection must be reproducable, specify a value for this argument.'
        },
    ],
    doc => 'Tool to merge BAM or SAM files aligned against a split reference (eg one that was larger than 4GiB before being broken up).'
};

sub help_detail {
    return 'Tool to merge BAM or SAM files aligned against a split reference (eg one that was larger than 4GiB before being broken up).  ' .
           'The alignments in the input files must be sorted by read name (samtools sort -n will do this).';
}

=pod

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
        # if softclip threshold is specified, we will include softclips in the number of mismatches
        my $rhs_mismatch = $rhs->[0]{NM} + $rhs->[1]{NM};
        my $lhs_mismatch = $lhs->[0]{NM} + $lhs->[1]{NM};
        if ($self->softclip_mismatch){
            $rhs_mismatch += $rhs->[0]{softclipped_bases} + $rhs->[1]{softclipped_bases}
            $lhs_mismatch += $lhs->[0]{softclipped_bases} + $lhs->[1]{softclipped_bases}
        }
        return $rhs_mismatch <=> $lhs_mismatch;
    }
    if($lhs_num_aln == 0) {
        # Neither alignment pair has a mapped read; the alignment pairs are equally bad print random pair
        return rand_select($rhs,$lhs);
    }

=cut

    
sub execute {
    my $self = shift;

    $self->dump_status_messages(1);
    $self->dump_error_messages(1);

    my @in_filenames = $self->input_files;
    my $out_filename = $self->output_file;

    if(scalar(@in_filenames) < 2) {
        die "Please supply at least two input files.\n";
    }

    srand($self->rand_seed);

    # Open input files and output file

    my @in_fh;
    foreach my $in_filename (@in_filenames) {
        my $in_fh = $self->open_bamsam_in($in_filename);
        die "failed to open file $in_filename! $!\n" unless $in_fh;
        push @in_fh, $in_fh;
    }

    my $out_fh = $self->open_bamsam_out($out_filename);

    my $in_count = @in_filenames;

    my $name_col     = 0;
    my $flag_col     = 1;
    my $ref_name_col = 2;
    my $ref_pos_col  = 3;
    my $map_qual_col = 4;
    my $cigar_col    = 5;

    my $name_count = 0;
    my $read_count = 0;
    my $finished_reading = 0;
    my $last_name = '';
    
    READ_NAME:
    for (;;) {
        # alignments for the next fragment of DNA, which may have one or two reads per file
        # $a[$file][$readnum][$col]
        my @a = ();      
        
        # the other metrics used to evaluate the alignment
        # currently there is just one: the NM tag holds edit distance, which is used to measure w
        # $m[$file][$readnum][$meric]
        my @m = ();

        my $name;
        my $paired_end1_count = 0;
        my $paired_end2_count = 0;
        my $single_end_count = 0;
        
        READ:
        for my $rn (1,2) {
            
            FILE:
            for (my $in_num=0; $in_num<$in_count; $in_num++) {
                my $fh = $in_fh[$in_num];
                my $line = <$fh>;

                if ( not defined $line ) {
                    $finished_reading++;
                    last READ_NAME;
                }

                # skip headers (pass through to output) 
                if(substr($line, 0, 1) eq '@') {
                    $out_fh->print($line);
                    redo;
                }

                # parse columns
                chomp $line;
                my @s = split(/\t/,$line);
                if(scalar(@s) < 10) {
                    die "Less than 10 columns in alignment record from $in_filenames[$in_num]: $line\n";
                }
                
                # ensure we are sorted by name
                if (not defined $name) {
                    $name = $s[$name_col];
                    if ($name eq $last_name) {
                        die "read $name appears in alignments after it is believed we are done processing that read's alignments!:\n$line";
                    }
                    # save this for the next read to ensure we are sorted and nothing odd is happening in the file
                    $last_name = $name;
                }
                elsif ($s[$name_col] ne $name) {
                    die "name mismatch, expected $name, got $s[$name_col] from file $in_filenames[$in_num]\n$line";
                }
                
                # If the alignment is mapped, put the NM tag it on the end of the array for the read
                # (we pop this off before printing)
                my $nm;
                unless($s[$flag_col] & 4) {
                    # mapped
                    if ($self->softclip_mismatch){
                        #if the softclip_mismatch arg is present, add the softclip(and hardclip) bases from the cigar string to the mismatch number
                        my $cigar = $s[$cigar_col];
                        my ($start_clip) = $cigar =~ /^(\d+)[hs]/i;
                        my ($stop_clip) = $cigar =~/(\d+)[hs]$/i;
                        $start_clip ||=0;
                        $stop_clip ||=0;
                        $nm += $start_clip+$stop_clip;
                    }
                    my $nm_pos = -1;
                    foreach my $tag (@s[11..$#s]) {
                        $nm_pos = index($tag,'NM:i:');
                        if ($nm_pos != -1) {
                            $nm += substr($tag,$nm_pos+5);
                            last;
                        }
                    }
                    if ($nm_pos == -1) {
                        die 'Alignment index ' . $s[$name_col] . ' in "' . $in_filenames[$in_num] . "\" is mapped but lacks an NM tag\n$line";
                    }
                }
                push @s, $nm; # this will be undef for unaligned reads

                # put the alignment in the buffer for this file & read position
                if ($s[$flag_col] & 64) {       # paired end 1
                    $a[$in_num][0] = \@s;
                    $paired_end1_count++;                    
                }
                elsif ($s[$flag_col] & 128 ) {  # paired end 2
                    $a[$in_num][1] = \@s;                   
                    $paired_end2_count++;
                }
                else {
                    $a[$in_num][0] = \@s;
                    $single_end_count++;
                }

            } # done checking each file for its next read

            $read_count++;

            last if $single_end_count;

        } # done getting each read for a given fragment


        #
        # process the read or pair across all files..        
        #


        # this is entirely to provide progress message.
        $name_count++;
        if ($name_count % 500_000 == 0) {
            print STDERR "fragments: $name_count, reads: $read_count\n";
        }        

        if ($finished_reading){
            print STDERR "Done reading! fragments: $name_count, reads: $read_count\n";
        }

        # sanity check
        if (
            ($paired_end1_count and $single_end_count)
                or
            ($paired_end1_count != $paired_end2_count)
            ) {
            die "Odd mix of fragment and paired data for read $name, p1,p2,se: $paired_end1_count, $paired_end2_count, $single_end_count";
            }

            # find the best alignment, or the list of best alignments if there is a tie among several
            my @best = ();

            # we may have a clear winner just by looking at number of reads mapped
            my $best_n_reads_mapped = 0;
            for (my $in_num=0; $in_num<$in_count; $in_num++) {
                my $candidate = $a[$in_num]; # examine the candidate read or pair from this file

                my $n_reads_mapped = 0;
                $n_reads_mapped++ if                        not $candidate->[0][$flag_col] & 4;
                $n_reads_mapped++ if $paired_end1_count and not $candidate->[1][$flag_col] & 4;
                if ($n_reads_mapped > $best_n_reads_mapped) {
                    @best = ($candidate);
                    $best_n_reads_mapped = $n_reads_mapped;
                }
                elsif ($n_reads_mapped == $best_n_reads_mapped) {
                    push @best, $candidate;
                }                
            }

            if ($best_n_reads_mapped == 0) {
                # no reads aligned, they are all equal, just dump the first pair 
                @best = ($best[0]);
            }

            # next prefer reads where both in a pair agree on reference
            if ($paired_end1_count and @best > 1) {
                my @matching_refs = 
                grep { $_->[0][$ref_name_col] ne '*' }               # skip pairs which agree on no reference
                grep { $_->[0][$ref_name_col] eq $_->[1][$ref_name_col] } 
                @best;

                if (@matching_refs and @matching_refs < @best) {
                    @best = @matching_refs;
                }
            }

            # next look for best edit distance (NM)  
            if (@best > 1) {
                # more than one is a top candidate looking just at alignment count
                # look at the edit distance
                my $best_nm_sum = 10000000; #pick a sufficiently high "worst" mismatch value
                for my $candidate (@best) {

                    my $nm_sum = 0;
                    $nm_sum += $candidate->[0][-1]    if                        not $candidate->[0][$flag_col] & 4;
                    $nm_sum += $candidate->[1][-1]    if $paired_end1_count and not $candidate->[1][$flag_col] & 4;                
                    if ($nm_sum < $best_nm_sum) {
                        @best = ($candidate);
                        $best_nm_sum = $nm_sum;
                    }
                    elsif ($nm_sum == $best_nm_sum) {
                        push @best, $candidate;
                    }                
                }            
            }

            # still more than one "best": select randomly
            if (@best > 1) {
                my $selected_n = int(rand(scalar(@best)));
                @best = ($best[$selected_n]);
            }

            # output the alignment or alignment pair
            for my $alignment (@{$best[0]}) {
                pop @$alignment;
                print $out_fh join("\t",@$alignment),"\n";    
            }

        } # done handling a given read name

        for my $fh ($out_fh,@in_fh) { $fh->close };
        return 1;
    }

    1;

