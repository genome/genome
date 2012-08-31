package Genome::Model::Tools::BioSamtools::CompareBams;

use strict;
use warnings;

use Genome;

class Genome::Model::Tools::BioSamtools::CompareBams {
    is  => 'Genome::Model::Tools::BioSamtools',
    has_input => [
        first_bam_file => {
            is => 'String',
            doc => 'First query sorted BAM file .',
        },
        second_bam_file => {
            is => 'String',
            doc => 'Second query sorted BAM file.',
        },
        merged_bam_file => {
            is => 'String',
            doc => 'The path to output the resulting merged, unsorted BAM file.',
        },
        alignment_stats_file => {
            is => 'String',
            doc => 'The path to output the stats file.',
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

    my $first_bam_file = $self->first_bam_file;
    my $second_bam_file = $self->second_bam_file;
    my $merged_bam_file = $self->merged_bam_file;
    my $output_fh = Genome::Sys->open_file_for_writing($self->alignment_stats_file);
    
    unless ($first_bam_file && $second_bam_file && $merged_bam_file) {
        die('Usage:  tophat_alignment_summary.pl <FASTQ_BAM> <TOPHAT_BAM> <MERGED_BAM>');
    }

    #my ($basename, $dirname, $suffix) = File::Basename::fileparse($merged_bam_file,qw/\.bam/);
    #my ($tmp_fh, $tmp_merged_bam_file) = tempfile($basename.'_XXXX',SUFFIX => $suffix, TMPDIR=>1);
    #$tmp_fh->close;

    my $merged_bam = Bio::DB::Bam->open($merged_bam_file,'w');
    unless ($merged_bam) {
        die('Failed to open output BAM file: '. $merged_bam_file);
    }

    my ($first_bam,$first_header) = validate_sort_order($first_bam_file);

    my ($second_bam,$second_header) = validate_sort_order($second_bam_file);
    &validate_second_bam_header($second_header);
    $merged_bam->header_write($second_header);
    
    my $target_names = $second_header->target_name;
    my %chr_hits;
    my $total_reads = 0;
    my $unmapped_count = 0;
    my $previous_second_read_name = '';
    while (my $second_read = $second_bam->read1) {
        my $second_flag = $second_read->flag;
        my $second_read_qname = $second_read->qname;
        my $second_read_end = 0;
        if ($second_flag & 1) {
            if ($second_flag & 64) {
                $second_read_end = 1;
            } elsif ($second_flag & 128) {
                $second_read_end = 2;
            } else {
                die ('Lost read pair info for: '. $second_read_qname);
            }
        }
        unless ($second_flag & 4) {
           # my $num_hits  = $second_read->aux_get('NH');
           # unless (defined($num_hits)) { die ('Failed to parse NH tag from BAM file: '. $second_bam_file); }
            # TODO: check mate chr and look for discordant pairs
            my $chr = $target_names->[$second_read->tid];
        my $second_read_name = $second_read_qname .'/'. $second_read_end;
        if ($second_read_name eq $previous_second_read_name) {
            next;
        }
        my $first_read = $first_bam->read1;
        my $first_flag = $first_read->flag;
        my $first_read_end;
        if ($first_flag & 64 ) {
            $first_read_end = 1;
        } elsif ( $first_flag & 128 ) {
            $first_read_end = 2;
        }
        my $first_read_name = $first_read->qname .'/'. $first_read_end;
        while ($first_read_name ne $second_read_name) {
            $unmapped_count++;
            $total_reads++;
            $merged_bam->write1($first_read);
            $first_read = $first_bam->read1;
            $first_flag = $first_read->flag;
            if ($first_flag & 64 ) {
                $first_read_end = 1;
            } elsif ( $first_flag & 128 ) {
                $first_read_end = 2;
            }
            $first_read_name = $first_read->qname .'/'. $first_read_end;
        }
        $total_reads++;
        $merged_bam->write1($second_read);
        $previous_second_read_name = $second_read_name;
    }

    # Write the remaining first reads
    while (my $first_read = $first_bam->read1) {
        $unmapped_count++;
        $total_reads++;
        $merged_bam->write1($first_read);
    }

    # sort by chr position putting unmapped reads at end of BAM
    # A call to samtools or picard may be faster/more efficient
    #Bio::DB::Bam->sort_core(0,$tmp_merged_bam_file,$dirname.'/'.$basename,4000);
    #unlink($tmp_merged_bam_file) || die('Failed to remove temp BAM file: '. $tmp_merged_bam_file);


    print $output_fh '##Total Reads: '. $total_reads ."\n";
    print $output_fh '##Unmapped Reads: '. $unmapped_count ."\n";
   # print $output_fh '##Unique Alignments: '. $total_top_hits ."\n";
 #   print $output_fh '##Total Reads Mapped: '. $total_mapped ."\n";
    $output_fh->close;
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

sub validate_second_bam_header {
    my $header = shift;
    my $text = $header->text;
    my @lines = split("\n",$text);
    my @pgs = grep { $_ =~ /^\@PG/ } @lines;
    unless (scalar(@pgs) == 1) {
        die('Found multiple PG lines in header: '. "\n\t". join("\n\t",@pgs) ."\nRefusing to continue parsing header.");
    }
    return 1;
}

}
1;
