package Genome::Model::Tools::BioSamtools::Tophat2AlignmentStats;

use strict;
use warnings;

use Genome;

class Genome::Model::Tools::BioSamtools::Tophat2AlignmentStats {
    is  => 'Genome::Model::Tools::BioSamtools',
    has_input => [
        bam_file => {
            is => 'String',
            doc => 'An aligned BAM file coordinate sorted with ALL alignments.',
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

    my $bam_file = $self->bam_file;
    my $output_fh = Genome::Sys->open_file_for_writing($self->alignment_stats_file);
    my $bam = Bio::DB::Bam->open($bam_file);
    unless ($bam) {
        die('Failed to open BAM file: '. $bam_file);
    }
    my $header = $bam->header;

    # The aligned BAM file should be output from Tophat
    $self->validate_aligned_bam_header($header);

    my $target_names = $header->target_name;

    my %chr_hits;
    my %total_hits;
    # Iterate over every alignment in the Tophat BAM file
    while (my $align = $bam->read1) {
        my $flag = $align->flag;
        my $chr = $target_names->[$align->tid];
        # Only count primary alignments
        if ($flag & 256) {
            if ($align->cigar_str =~ /N/) {
                $chr_hits{$chr}{'multi_spliced'}++;
                $total_hits{'multi_spliced'}++;
            }
            next;
        }

        $total_hits{reads}++;

        # Only work with mapped reads
        unless ($flag & 4) {
            $total_hits{mapped_count}++;
            # get the number of alignments per read
            my $num_hits  = $align->aux_get('NH');
            unless (defined($num_hits)) {
                die ('Failed to parse NH tag from BAM file: '. $bam_file);
            }

            # TODO: check mate chr and look for discordant pairs
            my $chr = $target_names->[$align->tid];

            # Uniquely aligned reads
            if ($num_hits == 1) {
                $chr_hits{$chr}{'top'}++;
                $total_hits{'top'}++;
                if ($align->cigar_str =~ /N/) {
                    $chr_hits{$chr}{'top_spliced'}++;
                    $total_hits{'top_spliced'}++;
                }
            # Multiple hits for read
            } elsif ($num_hits > 1) {
                $chr_hits{$chr}{'multi_reads'}++;
                $chr_hits{$chr}{'multi_hits'} += $num_hits;
                $total_hits{'multi_reads'}++;
                $total_hits{'multi_hits'} += $num_hits;
                if ($align->cigar_str =~ /N/) {
                    $chr_hits{$chr}{'multi_spliced'}++;
                    $total_hits{'multi_spliced'}++;
                }
            } else {
                die('No hits found for read '. $align->qname .' in bam file '. $bam_file);
            }
        } else {
            $total_hits{unmapped_count}++;
        }
    }

    print $output_fh "chr\ttop\ttop-spliced\tpct_top_spliced\tmulti_reads\tpct_multi_reads\tmulti_spliced_reads\tpct_multi_spliced_reads\tmulti_hits\n";
     for my $chr (sort keys %chr_hits) {
        my $top_hits = $chr_hits{$chr}{'top'} || 0;
        my $top_spliced = $chr_hits{$chr}{'top_spliced'} || 0;
        my $multi_reads = $chr_hits{$chr}{'multi_reads'} || 0;
        my $multi_hits = $chr_hits{$chr}{'multi_hits'} || 0;
        my $multi_spliced_reads = $chr_hits{$chr}{'multi_spliced'} || 0;

        my $pct_top_spliced = 0;
        if ($top_hits) {
            $pct_top_spliced = sprintf("%.02f",(($top_spliced / $top_hits) * 100));
        }

        my $pct_multi_hit_reads = 0;
        my $pct_multi_spliced_reads = 0;
        if ($multi_reads) {
            $pct_multi_hit_reads = sprintf("%.02f",( ( $multi_reads / ($multi_reads + $top_hits) ) * 100 ));
            $pct_multi_spliced_reads = sprintf("%.02f",( ( $multi_spliced_reads / ($multi_reads) ) * 100 ));
        }
        print $output_fh $chr ."\t". $top_hits ."\t". $top_spliced ."\t". $pct_top_spliced
            ."\t". $multi_reads ."\t". $pct_multi_hit_reads ."\t". $multi_spliced_reads
                ."\t". $pct_multi_spliced_reads ."\t". $multi_hits ."\n";
    }
    print $output_fh '##Total Reads: '. $total_hits{reads} ."\n";
    print $output_fh '##Unmapped Reads: '. $total_hits{unmapped_count} ."\n";
    print $output_fh '##Unique Alignments: '. $total_hits{top} ."\n";
    print $output_fh '##Multiple Hit Reads: '. $total_hits{'multi_reads'}."\n";
    print $output_fh '##Multiple Hit Sum: '. $total_hits{'multi_hits'}."\n";
    print $output_fh '##Total Reads Mapped: '. $total_hits{mapped_count} ."\n";
    $output_fh->close;

    my $expected_total_mapped = $total_hits{reads} - $total_hits{unmapped_count};
    unless ($expected_total_mapped == $total_hits{mapped_count}) {
        die('The expected number of mapped reads '. $expected_total_mapped .' does not match the sum '. $total_hits{mapped_count});
    }
    return 1;
}

sub validate_aligned_bam_header {
    my $self = shift;
    my $header = shift;
    
    my $text = $header->text;
    my @lines = split("\n",$text);
    my @pgs = grep { $_ =~ /^\@PG/ } @lines;
    my $found = 0;
    for my $pg_line (@pgs) {
        if ($pg_line =~/tophat/i) {
            $found = 1;
        }
    }
    unless ($found) {
        die('Input aligned BAM file '. $self->bam_file .' must be aligned with Tophat!');
    }
    return 1;
}

sub parse_alignment_stats_summary_hash_ref {
    my $class = shift;
    my $tophat_stats = shift;
    
    unless ($tophat_stats) {
        die('Must provide a tophat alignment stats file in '. __PACKAGE__);
    }
    
    my $tophat_fh = Genome::Sys->open_file_for_reading($tophat_stats);
    my %tophat_metrics;
    while (my $line = $tophat_fh->getline) {
        if ($line =~ /^##(.+):\s+(\d+)$/) {
            my $key = $1;
            my $value = $2;
            $key =~ s/ /_/g;
            $tophat_metrics{uc($key)} = $value;
        }
    }
    return \%tophat_metrics;
}

sub parse_alignment_stats_chromosome_hash_ref {
    my $class = shift;
    my $tophat_stats = shift;
    
    unless ($tophat_stats) {
        die('Must provide a tophat alignment stats file in '. __PACKAGE__);
    }
    my $tophat_fh = Genome::Sys->open_file_for_reading($tophat_stats);
    my ($tmp_fh,$tmp_file) = Genome::Sys->create_temp_file;
    while (my $line = $tophat_fh->getline) {
        if ($line !~ /^##(.+):\s+(\d+)$/) {
            print $tmp_fh $line;
        }
    }
    $tmp_fh->close;
    my $reader = Genome::Utility::IO::SeparatedValueReader->create(
        input =>$tmp_file,
        separator => "\t",
    );
    my %chr_metrics;
    while (my $data = $reader->next) {
        my $chr = delete($data->{chr});
        $chr_metrics{$chr} = $data;
    }
    return \%chr_metrics;
}

1;
