package Genome::Model::Tools::BioSamtools::TophatAlignmentStats;

use strict;
use warnings;

use Genome;

class Genome::Model::Tools::BioSamtools::TophatAlignmentStats {
    is  => 'Genome::Model::Tools::BioSamtools',
    has_input => [
        all_reads_bam_file => {
            is => 'String',
            doc => 'An unaligned BAM file querysorted with ALL reads from original FASTQ files.',
        },
        aligned_bam_file => {
            is => 'String',
            doc => 'A querysorted BAM file containing Tophat alignments.',
        },
        unaligned_bam_file => {
            is => 'String',
            doc => 'The path to output the resulting merged, unsorted BAM file.',
        },
        alignment_stats_file => {
            is => 'String',
            doc => 'A summary file of some calculated BAM alignment metrics.',
        },
    ],
};

sub help_synopsis {
    return <<EOS
    A Tophat based utility for alignment metrics.
EOS
}

sub help_brief {
    return <<EOS
    A Tophat based utility for alignment metrics.
EOS
}

sub help_detail {
    return <<EOS
--->Add longer docs here<---
EOS
}

sub execute {
    my $self = shift;

    my $all_reads_bam_file = $self->all_reads_bam_file;
    my $aligned_bam_file = $self->aligned_bam_file;
    my $unaligned_bam_file = $self->unaligned_bam_file;

    my $output_fh = Genome::Sys->open_file_for_writing($self->alignment_stats_file);

    my $unaligned_bam = Bio::DB::Bam->open($unaligned_bam_file,'w');
    unless ($unaligned_bam) {
        die('Failed to open output BAM file: '. $unaligned_bam_file);
    }

    # Bothe BAMs must be queryname sorted
    my ($all_reads_bam,$all_reads_header) = validate_sort_order($all_reads_bam_file);
    my ($aligned_bam,$aligned_header) = validate_sort_order($aligned_bam_file);

    # The aligned BAM file should be output from Tophat
    validate_aligned_bam_header($aligned_header);

    # The query names can be changed depending on the versions of software used.
    # If the query names do not match an expected format, fail until the exception is handled.
    unless (validate_query_name_format($all_reads_bam)) {
        die('Failed to validate query name format in unaligned BAM file: '. $all_reads_bam_file);
    }
    unless (validate_query_name_format($aligned_bam)) {
        die('Failed to validate query name format in aligned BAM file: '. $aligned_bam_file);
    }

    # Write the aligned BAM header to the output BAM
    $unaligned_bam->header_write($all_reads_header);

    my $target_names = $aligned_header->target_name;
    my %chr_hits;
    my %total_hits;
    my $previous_aligned_read_name = '';
    # Iterate over every read in the Tophat aligned BAM file
    while (my $aligned_read = $aligned_bam->read1) {
        my $aligned_flag = $aligned_read->flag;
        my $aligned_read_qname = $aligned_read->qname;
        # Determine the read end for paired-end sequence
        my $aligned_read_end = 0;
        if ($aligned_flag & 1) {
            if ($aligned_flag & 64) {
                $aligned_read_end = 1;
            } elsif ($aligned_flag & 128) {
                $aligned_read_end = 2;
            } else {
                die ('Lost read pair info for: '. $aligned_read_qname);
            }
        }

        # Add the read end to uniqify the query name;
        my $aligned_read_name = $aligned_read_qname .'/'. $aligned_read_end;

        # Only work with mapped reads
        unless ($aligned_flag & 4) {
            # get the number of alignments per read
            my $num_hits  = $aligned_read->aux_get('NH');
            unless (defined($num_hits)) { die ('Failed to parse NH tag from BAM file: '. $aligned_bam_file); }

            # TODO: check mate chr and look for discordant pairs
            my $chr = $target_names->[$aligned_read->tid];

            # Uniquely aligned reads
            if ($num_hits == 1) {
                $chr_hits{$chr}{'top'}++;
                $total_hits{'top'}++;
                if ($aligned_read->cigar_str =~ /N/) {
                    $chr_hits{$chr}{'top_spliced'}++;
                    $total_hits{'top_spliced'}++;
                }
            # Multiple hits for read
            } elsif ($num_hits > 1) {
                $chr_hits{$chr}{'multi'}{$aligned_read_name} = $num_hits;
                $total_hits{'multi'}{$aligned_read_name} = $num_hits;
                if ($aligned_read->cigar_str =~ /N/) {
                    $chr_hits{$chr}{'multi_spliced'}{$aligned_read_name}++;
                    $total_hits{'multi_spliced'}{$aligned_read_name}++;
                }
            } else {
                die('No hits found for '. $aligned_read_qname .'!  Please add support for unaligned reads.');
            }
        # This should never happen
        } else {
            die('Please add support for unaligned read found in Tophat BAM: '. $aligned_bam_file);
        }

        # This avoids duplicating an aligned read in the total read counts
        # Also there is no need to get an unaligned read until we hit a new aligned read name
        if ($aligned_read_name eq $previous_aligned_read_name) {
            next;
        }

        my $next_read = $all_reads_bam->read1;
        my $next_read_flag = $next_read->flag;
        my $next_read_end = 0;
        if ($next_read_flag & 1) {
            if ($next_read_flag & 64 ) {
                $next_read_end = 1;
            } elsif ( $next_read_flag & 128 ) {
                $next_read_end = 2;
            } else {
                die ('Lost read pair info for: '. $next_read->qname);
            }
        }
        my $next_read_name = $next_read->qname .'/'. $next_read_end;
        # If the unaligned read is not the same as the aligned read, we have an unmapped read
        # Loop over the unmapped reads until we get a name match, meaning it is aligned
        while ($next_read_name ne $aligned_read_name) {
            $total_hits{unmapped_count}++;
            $total_hits{reads}++;
            $unaligned_bam->write1($next_read);
            $next_read = $all_reads_bam->read1;
            $next_read_flag = $next_read->flag;
            if ($next_read_flag & 1) {
                if ($next_read_flag & 64 ) {
                    $next_read_end = 1;
                } elsif ( $next_read_flag & 128 ) {
                    $next_read_end = 2;
                } else {
                    die ('Lost read pair info for: '. $next_read->qname);
                }
            }
            $next_read_name = $next_read->qname .'/'. $next_read_end;
        }
        $total_hits{reads}++;
        $previous_aligned_read_name = $aligned_read_name;
    }

    # Write the remaining unaligned reads
    while (my $next_read = $all_reads_bam->read1) {
        $total_hits{unmapped_count}++;
        $total_hits{reads}++;
        $unaligned_bam->write1($next_read);
    }

    print $output_fh "chr\ttop\ttop-spliced\tpct_top_spliced\tmulti_reads\tpct_multi_reads\tmulti_spliced_reads\tpct_multi_spliced_reads\tmulti_hits\n";
    my $total_mapped;
    for my $chr (sort keys %chr_hits) {
        my $top_hits = $chr_hits{$chr}{'top'} || 0;
        $total_mapped += $top_hits;
        my $top_spliced = $chr_hits{$chr}{'top_spliced'} || 0;
        my $pct_top_spliced = 0;
        if ($top_hits) {
            $pct_top_spliced = sprintf("%.02f",(($top_spliced / $top_hits) * 100));
        }
        my $multi_reads = 0;
        my $multi_hits = 0;
        for my $read (keys %{$chr_hits{$chr}{'multi'}}) {
            $multi_reads++;
            $multi_hits += $chr_hits{$chr}{'multi'}{$read};
        }
        my $multi_spliced_reads = scalar(keys %{$chr_hits{$chr}{'multi_spliced'}}) || 0;
        my $pct_multi_hit_reads = 0;
        my $pct_multi_spliced_reads = 0;
        if ($multi_reads) {
            $pct_multi_hit_reads = sprintf("%.02f",( ( $multi_reads / ($multi_reads + $top_hits) ) * 100 ));
            $pct_multi_spliced_reads = sprintf("%.02f",( ( $multi_spliced_reads / ($multi_reads) ) * 100 ));
        }
        print $output_fh $chr ."\t". $top_hits ."\t". $top_spliced ."\t". $pct_top_spliced ."\t". $multi_reads ."\t". $pct_multi_hit_reads
            ."\t". $multi_spliced_reads ."\t". $pct_multi_spliced_reads ."\t". $multi_hits ."\n";
    }
    my $total_multi_reads = 0;
    my $total_multi_hits = 0;
    for my $read (keys %{$total_hits{'multi'}}) {
        $total_mapped++;
        $total_multi_reads++;
        $total_multi_hits += $total_hits{'multi'}{$read};
    }
    my $expected_total_mapped = $total_hits{reads} - $total_hits{unmapped_count};
    print $output_fh '##Total Reads: '. $total_hits{reads} ."\n";
    print $output_fh '##Unmapped Reads: '. $total_hits{unmapped_count} ."\n";
    print $output_fh '##Unique Alignments: '. $total_hits{top} ."\n";
    print $output_fh '##Multiple Hit Reads: '. $total_multi_reads ."\n";
    print $output_fh '##Multiple Hit Sum: '. $total_multi_hits ."\n";
    print $output_fh '##Total Reads Mapped: '. $total_mapped ."\n";
    $output_fh->close;
    unless ($expected_total_mapped == $total_mapped) {
        die('The expected number of mapped reads '. $expected_total_mapped .' does not match the sum '. $total_mapped);
    }
    return 1;
}


sub validate_sort_order {
    my $bam_file = shift;
    my $bam = Bio::DB::Bam->open($bam_file);
    unless ($bam) {
        die('Failed to open BAM file: '. $bam_file);
    }
    my $header = $bam->header;
    my $text = $header->text;
    my @lines = split("\n",$text);
    my @hds = grep { $_ =~ /^\@HD/ } @lines;
    unless (scalar(@hds) == 1) {
        die('Found multiple HD lines in header: '. "\n\t" . join("\n\t",@hds)) ."\nRefusing to continue parsing BAM file: ". $bam_file;
    }
    my $hd_line = $hds[0];
    if ($hd_line =~ /SO:(\S+)/) {
        my $sort_order = $1;
        unless ($sort_order eq 'queryname') {
            die('Input BAM files must be sorted by queryname!  BAM file found to be sorted by \''. $sort_order .'\' in BAM file: '. $bam_file);
        }
    } else {
        die('Input BAM files must be sorted by queryname!  No sort order found for input BAM file: '. $bam_file);
    }
    return ($bam,$header);
}

sub validate_aligned_bam_header {
    my $header = shift;
    my $text = $header->text;
    my @lines = split("\n",$text);
    my @pgs = grep { $_ =~ /^\@PG/ } @lines;
    unless (scalar(@pgs) == 1) {
        die('Found multiple PG lines in header: '. "\n\t". join("\n\t",@pgs) ."\nRefusing to continue parsing header.");
    }
    my $pg_line = $pgs[0];
    if ($pg_line =~ /ID:(\S+)/) {
        my $program = $1;
        unless ($program eq 'TopHat') {
            die('Input aligned BAM file must be aligned with Tophat!');
        }
    } else {
        die('Input aligned BAM file has no defined aligner program');
    }
    return 1;
}

sub validate_query_name_format {
    my $bam = shift;

    my $start_position = $bam->tell;
    my $align = $bam->read1;

    Bio::DB::Bam::seek($bam,$start_position,0);
    my $query_name = $align->qname;
    my $align_flag = $align->flag;
    if ($align_flag & 1) {
        # HWI-ST474_108856544:6:1:10000:107184#GCCNAT
        unless ($query_name =~ /^\S+:\d+:\d+:\d+:\d+[#ACTGN0]*$/) {
            warn('Query name '. $query_name .' is invalid!');
            return;
        }
    }
    return $query_name;
}

1;
