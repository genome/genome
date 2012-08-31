package Genome::Model::Tools::RefCov::Reference::Stats;

use strict;
use warnings;

use Genome;
use Statistics::Descriptive;

my $DEFAULT_WINDOW_SIZE = 10_000_000;

class Genome::Model::Tools::RefCov::Reference::Stats {
    has => [
        bam => {
            is => 'Bio::DB::Bam',
        },
        bam_index => {
            is => 'Bio::DB::Bam::Index',
        },
        window_size => {
            is_optional => 1,
            default_value => $DEFAULT_WINDOW_SIZE,
        },
        minimum_depth => {
            is_optional => 1,
            default_value => 0,
        },
    ],
    has_optional => {
        bed_file => {
            is => 'String',
            doc => 'A BED format file of regions to limit stats to',
        },
        _statistics => {},
    }
};

sub create {
    my $class = shift;
    my %params = @_;
    unless ($] > 5.010) {
        die "Bio::DB::Sam requires perl 5.10 or greater!";
    }
    require Bio::DB::Sam;
    # There seems to be a problem creating the UR object with these wrapped C classes
    my $bam = delete($params{bam});
    my $bam_index = delete($params{bam_index});
    my $self = $class->SUPER::create(%params);
    $self->bam($bam);
    $self->bam_index($bam_index);
    unless ($self->_load_stats) {
        die('Failed to load stats for BAM:  '. $self->bam);
    }
    return $self;
}

sub _load_stats {
    my $self = shift;

    my $bam = $self->bam;
    my $index = $self->bam_index;
    my $header = $bam->header();
    my $stats = Statistics::Descriptive::Sparse->new();
    my %target_name_index;
    my $i = 0;
    my $targets = $header->n_targets();
    my $target_names = $header->target_name();
    for my $target_name (@{ $target_names }) {
        $target_name_index{$target_name} = $i++;
    }# Make sure our index is not off
    unless ($targets == $i) {
        die 'Expected '. $targets .' targets but counted '. $i .' indices';
    }
    unless ($self->bed_file) {
        my $offset = $self->window_size;
        for (my $tid = 0; $tid < $targets; $tid++) {
            my $seq_id = $header->target_name->[$tid];
            my $length = $header->target_len->[$tid];
            for (my $start = 0; $start <= $length; $start += $offset) {
                my $end = $start + $offset;
                if ($end > $length) {
                    $end = $length;
                }
                my $coverage = $index->coverage($bam,$tid,$start,$end);
                if ($self->minimum_depth) {
                    my @coverage = @{$coverage};
                    my @sorted_coverage = sort {$a <=> $b} (@coverage);
                    my $depth = shift(@sorted_coverage);
                    while (defined($depth) && $depth < $self->minimum_depth) {
                        $depth = shift(@sorted_coverage);
                    }
                    if ($depth) {
                        unshift(@sorted_coverage,$depth);
                    }
                    $stats->add_data(@sorted_coverage);
                } else {
                    $stats->add_data($coverage);
                }
            }
        }
    } else {
        my $regions = Genome::Model::Tools::RefCov::ROI::Bed->create(
            file => $self->bed_file,
        );
        unless ($regions) {
            die('Failed to load BED region file '. $self->bed_file );
        }
        while (my $region = $regions->next_region) {
            my $id = $region->{name};
            my $target = $region->{chrom};
            my $start = $region->{start};
            my $end = $region->{end};
            my $length = $region->{length};

            # Here we get the $tid from the region chromosome
            my $tid = $target_name_index{$target};
            unless (defined $tid) { die('Failed to get tid for target '. $target); }
            
            # low-level API uses zero based coordinates
            # all regions should be zero based, but for some reason the correct length is never returned
            # the API must be expecting BED like inputs where the start is zero based and the end is 1-based
            # you can see in the docs for the low-level Bio::DB::BAM::Alignment class that start 'pos' is 0-based,but calend really returns 1-based
            my $coverage = $index->coverage( $bam, $tid, $start-1, $end);
            if ($self->minimum_depth) {
                my @coverage = @{$coverage};
                my @sorted_coverage = sort {$a <=> $b} (@coverage);
                my $depth = shift(@sorted_coverage);
                while (defined($depth) && $depth < $self->minimum_depth) {
                    $depth = shift(@sorted_coverage);
                }
                if ($depth) {
                    unshift(@sorted_coverage,$depth);
                }
                $stats->add_data(@sorted_coverage);
            } else {
                $stats->add_data($coverage);
            }
        }
    }
    $self->_statistics($stats);
    return 1;
}

sub mean_coverage {
    my $self = shift;
    my $stats = $self->_statistics;
    return $stats->mean;
}

sub standard_deviation{
    my $self = shift;
    my $stats = $self->_statistics;
    return $stats->standard_deviation;
}

sub bases_covered {
    my $self = shift;
    my $stats = $self->_statistics;
    return $stats->count;
}

sub variance {
    my $self = shift;
    my $stats = $self->_statistics;
    return $stats->variance;
}

sub minimum_coverage {
    my $self = shift;
    my $stats = $self->_statistics;
    return $stats->min;
}

sub maximum_coverage {
    my $self = shift;
    my $stats = $self->_statistics;
    return $stats->max;
}

1;
