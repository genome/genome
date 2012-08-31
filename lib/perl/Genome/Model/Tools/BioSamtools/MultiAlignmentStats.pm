package Genome::Model::Tools::BioSamtools::MultiAlignmentStats;

use strict;
use warnings;

use Genome;

my @DEFAULT_HEADERS = qw/
                            total_alignments
                            total_fragments
                            unmapped_fragments
                            pc_unmapped_fragments
                            mapped_fragments
                            pc_mapped_fragments
                            unique_alignment_fragments
                            pc_unique_alignment_fragments
                            multiple_alignment_fragments
                            pc_multiple_alignment_fragments
                            multiple_alignment_sum
                            multiple_alignment_per_fragment
                        /;

class Genome::Model::Tools::BioSamtools::MultiAlignmentStats {
    is  => 'Genome::Model::Tools::BioSamtools',
    has_input => [
        aligned_bam_file => {
            is => 'String',
            doc => 'A querysorted BAM file containing bwasw alignments.',
        },
        unique_bam_file => {
            is => 'String',
            is_optional => 1,
            doc => 'The path to output the resulting unique alignment BAM file.',
        },
        output_stats_tsv => {
            is => 'String',
            is_optional => 1,
            doc => 'The output stats tsv file name.  If not defined, stats print to STDOUT.',
        },
    ],
    has_optional => {
        _unique_bam => {},
        _aligned_bam => {},
        _total_alignments => { default_value => 0 },
        _unmapped_reads => { default_value => 0, },
        _unique_reads => { default_value => 0, },
        _multi_reads => { default_value => 0, },
        _multi_alignments => { default_value => 0, },
    },
};

sub help_synopsis {
    return <<EOS
    A bwasw based utility for alignment metrics.
EOS
}

sub help_brief {
    return <<EOS
    A bwasw based utility for alignment metrics.
EOS
}

sub help_detail {
    return <<EOS
Some aligners output multiple alignments per fragment in SAM/BAM format.
This tool is designed to take a queryname sorted, perferably by Picard, SAM/BAM file
and output metrics related to mapped/unmapped fragments and unique/multiple alignments.

Example of STDOUT:

Total Fragments: 940
Unmapped Fragments: 329	35.00%
Mapped Fragments: 611	65.00%
Unique Alignment Fragments: 300	31.91%
Multiple Hit Fragments: 311	33.09%
Multiple Hit Alignment Sum: 995
Multiple Hit Mean per Fragment: 3.20

Example tsv headers:

total_fragments
unmapped_fragments
pc_unmapped_fragments
mapped_fragments
pc_mapped_fragments
unique_alignment_fragments
pc_unique_alignment_fragments
multiple_alignment_fragments
pc_multiple_alignment_fragments
multiple_alignment_sum
multiple_alignment_per_fragment

EOS
}

sub execute {
    my $self = shift;

    # Alignment BAM must be queryname sorted
    my ($aligned_bam,$aligned_header) = $self->validate_sort_order('queryname');
    $self->_aligned_bam($aligned_bam);
    
    # Open unique BAM file if defined
    my $unique_bam_file = $self->unique_bam_file;
    my $unique_bam = undef;
    if ($unique_bam_file) {
        $unique_bam = Bio::DB::Bam->open($unique_bam_file,'w');
        unless ($unique_bam) {
            die('Failed to open output BAM file: '. $unique_bam_file);
        }
        $unique_bam->header_write($aligned_header);
        $self->_unique_bam($unique_bam);
    }

    $self->evaluate_by_read_name;
    
    # Print stats
    if ($self->output_stats_tsv) {
        $self->print_tsv;
    } else {
        $self->print_stdout;
    }
    return 1;
}

sub evaluate_by_read_name {
    my $self = shift;

    # Eval alignment objects
    my $previous_aligned_read = $self->read1;

    # Iterate over the remaining alignments
    while (my $aligned_read = $self->read1) {
        if ( $self->compare_reads($aligned_read,$previous_aligned_read) ) {
            # Multiple Alignments per Read
            my @alignments = ($previous_aligned_read);
            while ($aligned_read && $self->compare_reads($aligned_read,$previous_aligned_read) ) {
                push @alignments, $aligned_read;
                $previous_aligned_read = $aligned_read;
                $aligned_read = $self->read1;
            }
            $self->tally_alignments(@alignments);
        } else {
            # Unique Alignment for Read
            $self->tally_alignments($previous_aligned_read);
        }
        $previous_aligned_read = $aligned_read;
    }

    # Write last alignment if unique
    if ($previous_aligned_read) {
        $self->tally_alignments($previous_aligned_read);
    }
    
    return 1;
}

sub compare_reads {
    my $self = shift;
    my $query_read_a = shift;
    my $query_read_b = shift;

    if ($query_read_a->qname eq $query_read_b->qname) {
        my $query_read_a_flag = $query_read_a->flag;
        my $query_read_b_flag = $query_read_b->flag;
        if ( ($query_read_a_flag & 1) &&  ($query_read_b_flag & 1) ) {
            if ( ($query_read_a_flag & 64) && ($query_read_b_flag & 64) ) {
                # Both are read1
                return 1;
            } elsif ( ($query_read_a_flag & 128) && ($query_read_b_flag & 128) ) {
                # Both are read2
                return 1;
            }
        } else {
            # Same name for fragment read
            return 1;
        }
    }
    return 0;
}

sub validate_sort_order {

    my $self = shift;
    my $expected_sort_order = shift;
    my $bam_file = $self->aligned_bam_file;
    
    my $bam = Bio::DB::Bam->open($bam_file);
    unless ($bam) {
        $self->error_message('Failed to open BAM file: '. $bam_file);
        die($self->error_message);
    }
    my $header = $bam->header;
    my $text = $header->text;
    my @lines = split("\n",$text);
    my @hds = grep { $_ =~ /^\@HD/ } @lines;
    if (@hds) {
        unless (scalar(@hds) == 1) {
            $self->error_message('Found multiple HD lines in header: '. "\n\t" . join("\n\t",@hds) ."\nRefusing to continue parsing BAM file: ". $bam_file);
            die($self->error_message);
        }
        my $hd_line = $hds[0];
        if ($hd_line =~ /SO:(\S+)/) {
            my $sort_order = $1;
            unless ($sort_order eq $expected_sort_order) {
                $self->error_message('Input BAM files must be sorted by '. $expected_sort_order .'!  BAM file found to be sorted by \''. $sort_order .'\' in BAM file: '. $bam_file);
                die($self->error_message);
            }
        } else {
            $self->error_message('Input BAM files must be sorted by '. $expected_sort_order .'!  No sort order found for input BAM file: '. $bam_file);
            die($self->error_message);
        }
    } else {
        $self->warning_message('Failed to validate sort order!.  Assuming '. $expected_sort_order .' sorted!');
    }
    return ($bam,$header);
}

sub tally_alignments {
    my $self = shift;
    my @qname_alignments = @_;
    # This was almost refactored to parse read1 and read2 out but above we rely on the sorting
    my $top_alignment = $self->resolve_top_alignment(@qname_alignments);
    return $top_alignment;
}

sub resolve_top_alignment {
    my $self = shift;
    my @alignments = @_;

    unless (scalar(@alignments)) {
        die('Failed to resolve top alignment when no alignments passed to subroutine!');
    }
    
    my $top_alignment;
    for my $alignment (@alignments) {
        if ( !($alignment->flag & 256) ) {
            if (!defined($top_alignment)) {
                $top_alignment = $alignment;
            } elsif ($alignment->qual > $top_alignment->qual) {
                $top_alignment = $alignment;
            } elsif ($alignment->qual == $top_alignment->qual) {
                my $alignment_as_tag = $alignment->aux_get('AS');
                my $top_alignment_as_tag = $top_alignment->aux_get('AS');
                unless (defined($alignment_as_tag)) {
                    $self->warning_message('Missing AS tag for '. $alignment->qname);
                    next;
                }
                unless (defined($top_alignment_as_tag)) {
                    $self->warning_message('Missing AS tag for '. $top_alignment->qname);
                    next;
                }
                if ($alignment_as_tag > $top_alignment_as_tag) {
                    $top_alignment = $alignment;
                } elsif ($alignment_as_tag == $top_alignment_as_tag ) {
                    # TODO: determine best method for picking random alignment
                    # TODO: use seed to make rand consistent
                    $top_alignment = $alignments[rand @alignments];
                    # TODO: set mapping quality to zero ?
                }
            }
        }
    }
    if (scalar(@alignments) > 1) {
        $self->multi_mapped(scalar(@alignments));
    } else {
        $self->unique;
    }
    my $unique_bam = $self->_unique_bam;
    if ($unique_bam) {
        $unique_bam->write1($top_alignment);
    }
    return $top_alignment;
}

sub add_total_alignment {
    my $self = shift;
    my $total = $self->_total_alignments;
    $total++;
    $self->_total_alignments($total);
    unless ($total % 1_000_000) {
        $self->status_message('Processed '. $total .' alignments...');
    }
    return 1;
}

sub add_unmapped_read {
    my $self = shift;
    my $unmapped = $self->_unmapped_reads;
    $unmapped++;
    $self->_unmapped_reads($unmapped);
    return 1;
}

sub unique {
    my $self = shift;
    my $unique = $self->_unique_reads;
    $unique++;
    $self->_unique_reads($unique);
    return 1;
}

sub multi_mapped {
    my $self = shift;
    my $alignments = shift;

    my $multi = $self->_multi_reads;
    $multi++;
    $self->_multi_reads($multi);

    my $multi_align = $self->_multi_alignments;
    $multi_align += $alignments;
    $self->_multi_alignments($multi_align);
    return 1;
}

sub total_count {
    my $self = shift;
    my $total_reads = ($self->_unmapped_reads + $self->_unique_reads + $self->_multi_reads);
    return $total_reads;
}

sub mapped_count {
    my $self = shift;
    my $mapped_reads = ($self->_unique_reads + $self->_multi_reads);
    return $mapped_reads;
}

sub print_stdout {
    my $self = shift;

    my $total_reads = $self->total_count;
    my $mapped_reads = $self->mapped_count;
    print 'Total Fragments: '. $total_reads ."\n";
    print 'Unmapped Fragments: '. $self->_unmapped_reads ."\t". sprintf("%.02f",( ($self->_unmapped_reads / $total_reads) * 100)) ."%\n";
    print 'Mapped Fragments: '. $mapped_reads ."\t". sprintf("%.02f",( ($mapped_reads / $total_reads) * 100))  ."%\n";
    print 'Unique Alignment Fragments: '. $self->_unique_reads ."\t". sprintf("%.02f",( ($self->_unique_reads / $total_reads) * 100)) ."%\n";
    print 'Multiple Hit Fragments: '. $self->_multi_reads ."\t". sprintf("%.02f",( ($self->_multi_reads / $total_reads) * 100)) ."%\n";
    print 'Multiple Hit Alignment Sum: '.  $self->_multi_alignments ."\n";
    print 'Multiple Hit Mean per Fragment: '.  sprintf("%.02f",($self->_multi_alignments / $self->_multi_reads ))  ."\n";
    print 'Total Alignments: '. $self->_total_alignments ."\n";

    return 1;
}

sub headers {
    my $class = shift;
    return \@DEFAULT_HEADERS;
}

sub print_tsv {
    my $self = shift;
    my $tsv_writer = Genome::Utility::IO::SeparatedValueWriter->create(
        headers => $self->headers,
        separator => "\t",
        output => $self->output_stats_tsv,
    );
    unless ($tsv_writer) { die('Failed to open tsv writer: '. $self->output_stats_tsv) };
    my $total_reads = $self->total_count;
    my $mapped_reads = $self->mapped_count;
    my %data = (
        total_alignments => $self->_total_alignments,
        total_fragments => $total_reads,
        unmapped_fragments => $self->_unmapped_reads,
        pc_unmapped_fragments => sprintf("%.02f",( ($self->_unmapped_reads / $total_reads) * 100)),
        mapped_fragments => $mapped_reads,
        pc_mapped_fragments => sprintf("%.02f",( ($mapped_reads / $total_reads) * 100)),
        unique_alignment_fragments => $self->_unique_reads,
        pc_unique_alignment_fragments => sprintf("%.02f",( ($self->_unique_reads / $total_reads) * 100)),
        multiple_alignment_fragments => $self->_multi_reads,
        pc_multiple_alignment_fragments => sprintf("%.02f",( ($self->_multi_reads / $total_reads) * 100)),
        multiple_alignment_sum => $self->_multi_alignments,
        multiple_alignment_per_fragment =>  sprintf("%.02f",($self->_multi_alignments / $self->_multi_reads )),
    );
    $tsv_writer->write_one(\%data);
    $tsv_writer->output->close;
    return 1;
}

sub read1 {
    my $self = shift;
    my $aligned_read = $self->_aligned_bam->read1;
    $self->add_total_alignment();
    while ( $aligned_read && ($aligned_read->flag & 4) ) {
        $self->add_unmapped_read;
        $aligned_read = $self->read1;
    }
    return $aligned_read;
}


1;
