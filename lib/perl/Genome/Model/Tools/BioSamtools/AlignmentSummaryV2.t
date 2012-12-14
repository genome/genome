#!/usr/bin/env perl

use strict;
use warnings;

use Test::More;

use above 'Genome';
use Genome::Model::Tools::BioSamtools::AlignmentSummaryV2;

plan tests => 432;

# NOTE: When creating test data, remember that ROI start coords are 0-based,
# and the end coord is open-ended and 1-based.  For example start: 1000, end: 1001 is a region 1-long,
# and includes position 1000 but not 1001.
# BAM-file data coords are 0-based!
# a read's start position is the pos of the first mapped base (M or I), not any preceding S, H, P, etc

my $bam_reader;
my $bed_reader;
my $stats;

my @test_data = (
{
    # Test: 0
    # Read:    MMMM
    # ROI:           ----
    reads => [ 
        {   tid => 1,
            start => 1001,
            seq   => 'AAAA',
            cigar => '4M',
        } ],
    rois => [
        {  tid => 1,
            start => 1101,
            end => 1105,
         } ],
    stats => {
        total_bp => 4,
        unique_off_target_aligned_bp => 4,
        total_aligned_bp => 4,
        total_off_target_aligned_bp => 4,
    }
},

{
    # Test: 1
    # Read:          MMMM
    # ROI:    ----
    reads => [
        {   tid => 1,
            start => 1101,
            seq   => 'AAAA',
            cigar => '4M',
        } ],
    rois => [
        {  tid => 1,
            start => 1001,
            end => 1005,
        } ],
    stats => {
        total_bp => 4,
        unique_off_target_aligned_bp => 4,
        total_aligned_bp => 4,
        total_off_target_aligned_bp => 4,
    }
},
{
    # Test: 2
    # Read:  MMMM           MMMM
    # ROI:           ----
    reads => [
        {   tid => 1,
            start => 1001,
            seq   => 'AAAA',
            cigar => '4M',
        },
        {   tid => 1,
            start => 2001,
            seq   => 'AAAA',
            cigar => '4M',
        } ],
    rois => [
         {  tid => 1,
            start => 1101,
            end => 1105,
         } ],
    stats => {
        total_bp => 8,
        unique_off_target_aligned_bp => 8,
        total_aligned_bp => 8,
        total_off_target_aligned_bp => 8
    }
},

{
    # Test: 3
    # Read:          MMMM
    # ROI:  ----              ----
    reads => [
        {   tid => 1,
            start => 1101,
            seq   => 'AAAA',
            cigar => '4M',
        } ],
    rois => [
        {   tid => 1,
            start => 1001,
            end => 1005,
        },
        {   tid => 1,
            start => 2001,
            end => 1005,
        },
    ],
    stats => {
        total_bp => 4,
        unique_off_target_aligned_bp => 4,
        total_aligned_bp => 4,
        total_off_target_aligned_bp => 4,
    }
},
{
    # Test: 4
    # Read: MMMM
    # ROI:  ----  
    reads => [
        {   tid => 1,
            start => 1001,
            seq   => 'AAAA',
            cigar => '4M',
        } ],
    rois => [
         {  tid => 1,
            start => 1001,
            end => 1005,
         } ],
    stats => {
        total_bp => 4,
        unique_target_aligned_bp => 4,
        total_aligned_bp => 4,
        total_target_aligned_bp => 4,
    },
},
{
    # Test: 5
    # Read:   MMMM
    # ROI:  ----  
    reads => [
        {   tid => 1,
            start => 1000,
            seq   => 'AAAA',
            cigar => '4M',
        } ],
    rois => [
         {  tid => 1,
            start => 998,
            end => 1002,
         } ],
    stats => {
        total_bp => 4,
        unique_target_aligned_bp => 2,
        unique_off_target_aligned_bp => 2,
        total_aligned_bp => 4,
        total_target_aligned_bp => 2,
        total_off_target_aligned_bp => 2,
    }
},
{
    # Test: 6
    # Read: MMMM
    # ROI:    ----  
    reads => [
        {   tid => 1,
            start => 1000,
            seq   => 'AAAA',
            cigar => '4M',
        } ],
    rois => [
         {  tid => 1,
            start => 1002,
            end => 1006,
         } ],
    stats => {
        total_bp => 4,
        unique_target_aligned_bp => 2,
        unique_off_target_aligned_bp => 2,
        total_aligned_bp => 4,
        total_target_aligned_bp => 2,
        total_off_target_aligned_bp => 2,
    }
},

{
    # Test: 7
    # Read:    MMMM
    # ROI:  -----------
    reads => [
        {   tid => 1,
            start => 1000,
            seq   => 'AAAA',
            cigar => '4M',
        } ],
    rois => [
         {  tid => 1,
            start => 997,
            end => 1007,
         } ],
    stats => {
        total_bp => 4,
        unique_target_aligned_bp => 4,
        total_aligned_bp => 4,
        total_target_aligned_bp => 4,
    }
},

{
    # Test: 8
    # Read:   MMMM
    # Roi:     -
    reads => [
        {   tid => 1,
            start => 1000,
            seq   => 'AAAA',
            cigar => '4M',
        } ],
    rois => [
         {  tid => 1,
            start => 1001,
            end => 1002,
         } ],
    stats => {
        total_bp => 4,
        unique_target_aligned_bp => 1,
        total_aligned_bp => 4,
        total_target_aligned_bp => 1,
        unique_off_target_aligned_bp => 3,
        total_off_target_aligned_bp => 3,
    } 
},
{
    # Test: 9
    # read:   SSSSMMMM
    # roi:      ----
    reads => [
        {   tid => 1,
            start => 1004,
            seq   => 'AAAAAAAA',
            cigar => '4S4M',
        } ],
    rois => [
        {   tid => 1,
            start => 1002,
            end => 1006,
        } ],
    stats => {
        total_bp => 8,
        unique_target_aligned_bp => 2,
        total_aligned_bp => 4,
        total_target_aligned_bp => 2,
        unique_off_target_aligned_bp => 2,
        total_off_target_aligned_bp => 2,
        total_soft_clipped_bp => 4,
        total_unaligned_bp => 4,
    }
},

{
    # Test: 10
    # Read:   MMMMMMMMMSSSSSSSSSS
    # Roi:                -- --
    reads => [
        {   tid => 1,
            start => 1000,
            seq   => 'A'x100,
            cigar => '87M13S',
        } ],
    rois => [
        {   tid => 1,
            start => 1090,
            end => 1092,
        },
        {   tid => 1,
            start => 1095,
            end => 1097,
        } ],
    stats => {
        total_bp => 100,
        unique_target_aligned_bp => 0,
        total_aligned_bp => 87,
        total_target_aligned_bp => 0,
        unique_off_target_aligned_bp => 87,
        total_off_target_aligned_bp => 87,
        total_soft_clipped_bp => 13,
        total_unaligned_bp => 13,
    }
},

{
    # Test: 11
    # Read:     MMMMMMMMMMMMMMMMMMMM
    # Roi:    --------      ------------
    reads => [
        {   tid => 1,
            start => 1000,
            seq   => 'A'x20,
            cigar => '20M',
        } ],
    rois => [
        {   tid => 1,
            start => 998,
            end => 1006,
        },
        {   tid => 1, 
            start => 1012,
            end => 1024, 
        } ],
    stats => {
        total_bp => 20,
        unique_target_aligned_bp => 14,
        total_aligned_bp => 20,
        total_target_aligned_bp => 14,
        unique_off_target_aligned_bp => 6,
        total_off_target_aligned_bp => 6,
    }
},
{
    # Test: 12
    # Read1:     MMMM
    # Read2:     MMMM
    # Roi:           ----
    reads => [
        {   tid => 1,
            start => 1000,
            seq   => 'AAAA',
            cigar => '4M',
        },
        {   tid => 1,
            start => 1000,
            seq   => 'AAAA',
            cigar => '4M',
            is_duplicate => 1,
        } ],
    rois => [
        {   tid => 1,
            start => 1005,
            end => 1009,
        },],
    stats => {
        total_bp => 8,
        total_aligned_bp => 8,
        total_duplicate_bp => 4,
        unique_off_target_aligned_bp => 4,
        duplicate_off_target_aligned_bp => 4,
        total_off_target_aligned_bp => 8,
    }
},
{
    # Test: 13
    # Read:    MMMM
    # ROI:           ----
    reads => [ 
        {   tid => 1,
            start => 1001,
            seq   => 'AAAA',
            cigar => '4M',
            is_unmapped => 1
        } ],
    rois => [
        {  tid => 1,
            start => 1101,
            end => 1105,
         } ],
    stats => {
        total_bp => 4,
        total_unaligned_bp => 4,
    }
},
{
    # Test: 14
    # Read1:     MMMM
    # Read2:     MMMM
    # Roi:       ----
    reads => [
        {   tid => 1,
            start => 1000,
            seq   => 'AAAA',
            cigar => '4M',
        },
        {   tid => 1,
            start => 1000,
            seq   => 'AAAA',
            cigar => '4M',
            is_duplicate => 1,
        } ],
    rois => [
        {   tid => 1,
            start => 1000,
            end => 1005,
        },],
    stats => {
        total_bp => 8,
        total_aligned_bp => 8,
        total_duplicate_bp => 4,
        unique_target_aligned_bp => 4,
        duplicate_target_aligned_bp => 4,
        total_target_aligned_bp => 8,
    }
},
{
    # Test 15
    # Read:     MMMMIMMM
    # roi:    ------------
    reads => [
        {   tid => 1,
            start => 1000,
            seq   => 'AAAAAAAA',
            cigar => '4M1I3M',
        } ],
    rois => [
        {   tid => 1,
            start => 900,
            end => 1100,
        } ],
    stats => {
        total_bp => 8,
        total_aligned_bp => 8,
        unique_target_aligned_bp => 8,
        total_target_aligned_bp => 8,
    }
},
{
    # Test 16
    # Read:          SSSSMMMM
    # roi:         ----
    reads => [
        {   tid => 1,
            start => 1004,
            seq   => 'AAAAAAAA',
            cigar => '4S4M',
        } ],
    rois => [
        {   tid => 1,
            start => 998,
            end => 1002,
        } ],
    stats => {
        total_bp => 8,
        total_aligned_bp => 4,
        unique_target_aligned_bp => 0,
        total_target_aligned_bp => 0,
        total_unaligned_bp => 4,
        total_off_target_aligned_bp => 4,
        total_soft_clipped_bp => 4,
        unique_off_target_aligned_bp => 4,
    }
},
{
    # Test 17
    # Read:      MMMM
    # ROI:   ----
    reads => [
        {   tid => 1,
            start => 1000,
            seq   => 'AAAA',
            cigar => '4M',
        }],
    rois => [
        {   tid => 1,
            start => 996,
            end  => 1000,
        }],
    stats => {
        total_bp => 4,
        total_aligned_bp => 4,
        unique_off_target_aligned_bp => 4,
        total_off_target_aligned_bp => 4,
    },
},
{
    # test 18
    # Read:    MMMIMMMIMMM
    # ROI:  --------------
    reads => [
        {   tid => 1,
            start => 1000,
            seq   => 'AAACAAACAAA',
            cigar => '3M1I3M1I3M',
        }],
    rois => [
        {   tid => 1,
            start => 997,
            end   => 1011,
        }],
    stats => {
        total_bp => 11,
        total_aligned_bp => 11,
        unique_target_aligned_bp => 11,
        total_target_aligned_bp  => 11,
    },
},
{
    # test 19
    # Read: SSSSMMMM
    # Roi:    ----
    reads => [
        {   tid => 1,
            start => 1000,  # Note that this is the position of the first M, not the S
            seq => 'AAAAAAAA',
            cigar => '4S4M',
        }],
    rois => [
        {   tid => 1,
            start => 998,
            end => 1002,
        }],
    stats => {
        total_bp => 8,
        total_aligned_bp => 4,
        unique_target_aligned_bp => 2,
        total_target_aligned_bp => 2,
        unique_off_target_aligned_bp => 2,
        total_off_target_aligned_bp => 2,
        total_unaligned_bp => 4,
        total_soft_clipped_bp => 4,
    },
},
{
    # test 20
    # read: SSSSMMMM
    # Roi:      ----
    reads => [
        {   tid => 1,
            start => 1000,  # Note that this is the position of the first M, not the S
            seq => 'AAAAAAAA',
            cigar => '4S4M',
        }],
    rois => [
        {   tid => 1,
            start => 1000,
            end => 1004,
        }],
    stats => {
        total_bp => 8,
        total_aligned_bp => 4,
        unique_target_aligned_bp => 4,
        total_target_aligned_bp => 4,
        total_unaligned_bp => 4,
        total_soft_clipped_bp => 4,
    },
},
{
    # test 21
    # read: SSSSMMMM
    # Roi:        ----
    reads => [
        {   tid => 1,
            start => 1000,  # Note that this is the position of the first M, not the S
            seq => 'AAAAAAAA',
            cigar => '4S4M',
        }],
    rois => [
        {   tid => 1,
            start => 1002,
            end => 1006,
        }],
    stats => {
        total_bp => 8,
        total_aligned_bp => 4,
        unique_target_aligned_bp => 2,
        total_target_aligned_bp => 2,
        unique_off_target_aligned_bp => 2,
        total_off_target_aligned_bp => 2,
        total_unaligned_bp => 4,
        total_soft_clipped_bp => 4,
    },
},
{
    # test 22
    # read: SSSSMMMM
    # roi:    --------
    reads => [
        {   tid => 1,
            start => 1000,  # Note that this is the position of the first M, not the S
            seq => 'AAAAAAAA',
            cigar => '4S4M',
        }],
    rois => [
        {   tid => 1,
            start => 998,
            end => 1006,
        }],
    stats => {
        total_bp => 8,
        total_aligned_bp => 4,
        unique_target_aligned_bp => 4,
        total_target_aligned_bp => 4,
        total_unaligned_bp => 4,
        total_soft_clipped_bp => 4,
    },
},
{
    reads => [
        {   tid => 1,
            start => 1000,
            seq => 'AAAAA',
            cigar => '3M2N2M',
        }],
    rois => [
        {   tid => 1,
            start => 1000,
            end => 1006,
        }],
    stats => {
        total_bp => 5,
        total_aligned_bp => 5,
        unique_target_aligned_bp => 5,
        total_target_aligned_bp => 5,
    },
},
);

for (my $i = 0; $i < @test_data; $i++) {
    my $bam_reader = FakeBamReader->create(@{$test_data[$i]->{reads}});
    my $bed_reader = create_bed_reader(@{$test_data[$i]->{rois}});
    my $stats = Genome::Model::Tools::BioSamtools::AlignmentSummaryV2->summarize_alignments($bam_reader, $bed_reader);
    ok($stats, "Get alignment summary stats for test $i");
    foreach my $key ( keys %{$test_data[$i]->{stats}} ) {
        my $expected_value = $test_data[$i]->{stats}->{$key} || 0;
        my $value = delete($stats->{$key}) || 0;
        is($value, $expected_value, "$key is $value");
    }
    is($stats->{$_} || 0, 0, "$_ is 0") foreach keys %$stats;
}




sub create_bed_reader {
    my %rois;
    foreach my $roi ( @_ ) {
        my $tid = $roi->{'tid'};
        $rois{$tid} ||= [];
        push @{$rois{$tid}}, $roi;
    }

    foreach my $tid ( keys %rois ) {
        $rois{$tid} = [ sort { $a->{'start'} <=> $b->{'start'} } @{ $rois{$tid} } ];
    }

    return sub {
        my $tid = shift;
        my $roi = shift @{ $rois{$tid} };
        return unless $roi;
        return ( $roi->{'start'}, $roi->{'end'});
    }
}


# This package emulates Bio::DB::Sam objects
package FakeBamReader;

sub create {
    my $class = shift;

    my @reads = sort { $a->{'tid'} <=> $b->{'tid'} }
                sort { $a->{'start'} <=> $b->{'start'} }
                @_;

    my $self = [ map { FakeBamReader::Alignment->create($_) } @reads ];
    
    bless $self, $class;
    return $self;
}

sub read1 {
    my $self = shift;
    my $read = shift @$self;
    return unless $read;
    return $read->_clone;
}

package FakeBamReader::Alignment;
use Bio::DB::Sam::Constants qw(BAM_CIGAR_MASK BAM_CIGAR_SHIFT);

use constant cigar_ops => {'M' => 0, 'I' => 1, 'D' => 2, 'N' => 3, 'S' => 4, 'H' => 5, 'P' => 6, '=' => 7, 'X' => 8};
sub create {
    my $class = shift;
    my $read = shift;

    my $flag = 0;
    $flag |= 1024 if $read->{'is_duplicate'};
    $flag |= 4    if $read->{'is_unmapped'};
    $flag |= 1    if $read->{'is_paired'};
    $flag |= 2    if $read->{'is_proper_pair'};
    $flag |= 64   if $read->{'is_read1'};
    $flag |= 8    if $read->{'is_mate_unpapped'};
    $read->{'flag'} = $flag;

    my $cigar = [];
    my $cigar_string = delete $read->{'cigar'};
    $read->{'cigar'} = $cigar;
    while ($cigar_string =~ s/^(\d+)(\D+)//) {
        my $cigar_op = cigar_ops->{$2} & BAM_CIGAR_MASK;
        my $cigar_len = $1;
        my $data = $cigar_op | ( $cigar_len << BAM_CIGAR_SHIFT );
        push @$cigar, $data;
    }

    bless $read, $class;
    return $read;
}

sub flag {
    return shift->{'flag'};
}

sub l_qseq {
    return length(shift->{'seq'});
}

sub tid {
    return shift->{'tid'};
}

sub pos {
    return shift->{'start'};
}

sub calend {
    my $read = shift;
    return $read->{'start'} + length($read->{'seq'});
}

sub cigar {
    return [ @{ shift->{'cigar'} } ];
}

sub _clone {
    my $self = shift;
    my $copy = { %$self };
    bless $copy, ref $self;
    return $copy;
}
