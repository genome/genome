package Genome::Model::Tools::BioSamtools::ClusterCoverage;

use strict;
use warnings;

use Genome;
use PDL;
use PDL::NiceSlice;

my $DEFAULT_MINIMUM_ZENITH = '1';
my $DEFAULT_MINIMUM_DEPTH = '1';
my $DEFAULT_OFFSET = 50_000_000;
my $DEFAULT_MAXIMUM_DEPTH = 100_000_000;
my @DEFAULT_STATS_HEADERS = qw/
                                  name
                                  mean
                                  prms
                                  med
                                  min
                                  max
                                  adev
                                  rms
                              /;

class Genome::Model::Tools::BioSamtools::ClusterCoverage {
    is => 'Genome::Model::Tools::BioSamtools',
    has_input => [
        bam_file => {
            is => 'Text',
            doc => 'A path to a BAM format file of aligned capture reads',
        },
        minimum_zenith => {
            is => 'Text',
            doc => 'A minimum zenith(maximum depth) size to retain a cluster',
            default_value => $DEFAULT_MINIMUM_ZENITH,
            is_optional => 1,
        },
        minimum_depth => {
            is => 'Text',
            doc => 'A comma separated list of minimum depths to filter coverage',
            default_value => $DEFAULT_MINIMUM_DEPTH,
            is_optional => 1,
        },
        maximum_depth => {
            is => 'Text',
            doc => 'The maximum depth at one position for samtools to iterate through a pileup.',
            default_value => $DEFAULT_MAXIMUM_DEPTH,
            is_optional => 1,
        },
        offset => {
            is => 'Text',
            doc => 'The size(bp) of the sliding window offset.',
            default_value => $DEFAULT_OFFSET,
            is_optional => 1,
        },
        minimum_base_quality => {
            is => 'Text',
            doc => 'A minimum base quality to consider in coverage assesment.  THIS IS DEPRECATED FOR NOW.',
            is_deprecated => 1,
            is_optional => 1,
        },
        minimum_mapping_quality => {
            is => 'Text',
            doc => 'A minimum mapping quality to consider in coverage assesment.  THIS IS DEPRECATED FOR NOW.',
            is_deprecated => 1,
            is_optional => 1,
        },
        bed_file => {
            is => 'Text',
            doc => 'The output BED format file of clusters.',
        },
        stats_file => {
            is => 'Text',
            doc => 'Calculate statistics across clusters and print to file.',
            is_optional => 1,
        },
        maximum_depth_position => {
            is => 'Boolean',
            doc => 'Print the coordinates of the maximum depth positions in the stats output file.',
            is_optional => 1,
            default_value => 0,
        },
    ],
};

sub execute {
    my $self = shift;

    my $bed_fh = Genome::Sys->open_file_for_writing($self->bed_file);
    unless ($bed_fh) {
        die('Failed to open file for writing: '. $self->bed_file);
    }
    my $stats_writer;
    if ($self->stats_file) {
        my $headers = $self->resolve_stats_headers;
        $stats_writer = Genome::Utility::IO::SeparatedValueWriter->create(
            output => $self->stats_file,
            separator => "\t",
            headers => $headers,
        );
        unless ($stats_writer) {
            die('Failed to open file for writing: '. $self->stats_file);
        }
    }
    my $refcov_bam  = Genome::Model::Tools::RefCov::Bam->create(bam_file => $self->bam_file );
    unless ($refcov_bam) {
        die('Failed to load bam file '. $self->bam_file);
    }
    my $bam  = $refcov_bam->bio_db_bam;
    my $index = $refcov_bam->bio_db_index;
    my $header = $bam->header();

    # Number of reference sequences
    my $targets = $header->n_targets();

    # The reference sequence names in an array ref with indexed positions
    my $target_names = $header->target_name();

    # at the low level API the seq_id/target_name is meaningless
    # cache the target_names in a hash by actual reference sequence name
    # then we can look up the target index on the fly
    my %target_name_index;
    my $i = 0;
    for my $target_name (@{ $target_names }) {
        $target_name_index{$target_name} = $i++;
    }

    # Make sure our index is not off
    unless ($targets == $i) {
        die 'Expected '. $targets .' targets but counted '. $i .' indices';
    }

    my $quality_coverage_callback = sub {
        my ($tid,$pos,$pileups,$data) = @_;
        die('I have not implemented or tested filtering by mapping or base quality');
        my ($start,$end,$coverage) = @$data;
        #Here the position $pos is always zero-based, but the end position has to be 1-based in the coverage function
        if ($pos < $start || $pos >= $end) { return; }
        my $index = $pos - $start;
        for my $pileup (@$pileups) {
            my $base_position = $pileup->qpos;
            my $alignment = $pileup->alignment;
            if (defined($self->minimum_mapping_quality)) {
                unless ($alignment->qual >= $self->minimum_mapping_quality) {
                    next;
                }
            }
            my @base_qualities = $alignment->qscore;
            my $quality = $base_qualities[$base_position];
            if ($quality >= $self->minimum_base_quality) {
                @$coverage[$index]++;
            }
        }
    };
    my $offset = $self->offset;
    my $min_depth = $self->minimum_depth;
    for (my $tid = 0; $tid < $targets; $tid++) {
        my @previous_clusters;
        my $chr = $header->target_name->[$tid];
        #unless ($chr eq 'MT') { next; }
        #print 'Starting: '. $chr ."\n";
        my $chr_length = $header->target_len->[$tid];
        for (my $start = 0; $start <= $chr_length; $start += $offset) {
            my $end = $start + $offset;
            if ($end > $chr_length) {
                $end = $chr_length;
            }

            # low-level API uses zero based coordinates
            # all regions should be zero based, but for some reason the correct length is never returned
            # the API must be expecting BED like inputs where the start is zero based and the end is 1-based
            # you can see in the docs for the low-level Bio::DB::BAM::Alignment class that start 'pos' is 0-based,but calend really returns 1-based
            my $coverage;
            if (defined($self->minimum_base_quality) || defined($self->minimum_mapping_quality)) {
                #Start with an empty array of zeros
                my @coverage = map { 0 } (1 .. $offset);
                $coverage = \@coverage;
                # the pileup callback will add each base gt or eq to the quality_filter to the index position in the array ref
                $index->pileup($bam,$tid,$start,$end,$quality_coverage_callback,[$start,$end,$coverage])
            } else {
                # The bin size is set to zero to return each position
                $coverage = $index->coverage( $bam, $tid, $start, $end,0,$self->maximum_depth);
                #print $start ."\t". $end ."\n";
            }
            #my @clusters = $self->cluster_coverage($coverage,$min_depth,$start);
            #print 'Starting: '. $chr ."\t". $start ."\t". scalar(@{$coverage}) ."\n";
            my @clusters = $self->pdl_clusters($coverage,$min_depth,$start);
            if (scalar(@clusters)) {
                if (scalar(@previous_clusters)) {
                    my $last_cluster = $previous_clusters[-1];
                    my $first_cluster = $clusters[0];
                    #print Data::Dumper::Dumper($last_cluster, $first_cluster);
                    if (($first_cluster->[0] == $last_cluster->[1])) {
                        #Blunt end clusters: merge
                        $last_cluster->[1] = $first_cluster->[1];
                        #Merge the stats
                        if ($stats_writer) {
                            my $last_pdl = $last_cluster->[2];
                            my $first_pdl = $first_cluster->[2];
                            $last_pdl->append($first_pdl);
                        }
                        shift(@clusters);
                    }
                    $self->print_clusters($chr,\@previous_clusters,$bed_fh,$stats_writer);
                }
                @previous_clusters = @clusters;
            }
            #print 'Finished: '. $chr ."\t". $start ."\t". scalar(@{$coverage}) ."\n";
        }
        #print 'Printing remaining clusters: '. scalar(@previous_clusters) ."\n";
        $self->print_clusters($chr,\@previous_clusters,$bed_fh,$stats_writer);
        #print 'Finished: '. $chr ."\n";
    }
    $bed_fh->close;
    return 1;
}


sub pdl_clusters {
    my $self = shift;

    my $coverage = shift;
    my $min_depth = shift;
    my $offset = shift;

    unless (scalar(@{$coverage})) { return };
    #print Data::Dumper::Dumper($coverage) ."\n";
    my $chr_pdl = pdl $coverage;
    my ($quantity,$value) = rle($chr_pdl >= $min_depth);
    #print 'Quantity: '. $quantity->info ."\n";
    #print 'Value: '. $value->info ."\n";
    my ($padded_quantity,$padded_value) = $quantity->where($value,$quantity!=0);
    #print 'Padded Quantity: '. $padded_quantity ."\n";
    #print 'Padded Value: '. $padded_value ."\n";
    my $padded_quantity_cumsum = $padded_quantity->cumusumover;
    #print 'Padded Quantity cumsum: '. $padded_quantity_cumsum ."\n";
    my ($chr_gt_min_depth_idx,$chr_lt_min_depth_idx) = which_both($padded_value);
    # NO CLUSTERS
    unless ($chr_gt_min_depth_idx->nelem) {
        return;
    }

    my $start;
    my $end;
    # NO GAPS - CONTIGUOUS CLUSTER
    unless ($chr_lt_min_depth_idx->nelem) {
        # The start position is index zero
        $start = zeroes(1);
        # The end position is the cumsum
        $end = $padded_quantity_cumsum->index($chr_gt_min_depth_idx);
        #print 'Start: '. $start ."\n";
        #print 'End: '. $end ."\n";
    } else {
        my $gt_first = $chr_gt_min_depth_idx(0)->sclr;
        #print 'Gt_first: '. $gt_first ."\n";
        my $lt_first = $chr_lt_min_depth_idx(0)->sclr;
        #print 'Lt_first: '. $gt_first ."\n";

        my $gt_last = $chr_gt_min_depth_idx(-1)->sclr;
        #print 'Gt_last: '. $gt_last ."\n";
        my $lt_last = $chr_lt_min_depth_idx(-1)->sclr;
        #print 'Lt_last: '. $lt_last ."\n";
        $start = $padded_quantity_cumsum->index($chr_lt_min_depth_idx);
        $end = $padded_quantity_cumsum->index($chr_gt_min_depth_idx);
        if ($gt_first < $lt_first) {
            # Coverage first because the first index is zero
            # Add the initial start position(ie. zero)
            $start = append(zeroes(1),$start);
        }
        if ($gt_last < $lt_last) {
            # Gap on the end, remove the last value
            $start = $start(0:-2);
        }
    }
    my @clusters;
    # Iterate over each start position
    for (my $i = 0; $i < $start->getdim(0); $i++) {
        # Start is zero-based and End is one-based, just like BED format
        my $start_coordinate = $start($i)->sclr;
        # The end position is really the start of the gap
        my $end_coordinate = $end($i)->sclr - 1;
        my $cluster_pdl = $chr_pdl($start_coordinate:$end_coordinate)->sever;
        my ($gt_idx,$zero_idx) = which_both($cluster_pdl >= $min_depth);
        # TODO: This can be removed once everything is validated
        if ($zero_idx->getdim(0)) {
            die('Unexpected number in PDL: '. $start_coordinate ."\t". $end_coordinate ."\t". $cluster_pdl);
        }
        my $start_pos = $start_coordinate + $offset;
        my $end_pos = $end_coordinate + $offset + 1;
        push @clusters, [ $start_pos, $end_pos, $cluster_pdl];
    }
    return @clusters;
}

sub _round {
    my $self = shift;
    my $value = shift;
    return sprintf( "%.2f", $value );
}

sub print_clusters {
    my $self = shift;
    my $chr = shift;
    my $clusters = shift;
    my $bed_fh = shift;
    my $stats_writer = shift;

    for my $cluster (@{$clusters}) {
        my %data;
        my $start = $cluster->[0];
        my $end = $cluster->[1];
        my $pdl = $cluster->[2];
        my $name = $chr .':'. $start .'-'. $end;
        unless (defined($pdl)) {
            die('A cluster is defined with no coverage PDL object: '. Data::Dumper::Dumper($cluster));
        }
        my ($mean,$prms,$med,$min,$max,$adev,$rms) = $pdl->stats;

        if ($max >= $self->minimum_zenith) {
            print $bed_fh $chr ."\t". $start ."\t". $end ."\t". $name ."\n";
            if ($stats_writer) {
                %data = (
                    name => $name,
                    mean => $self->_round($mean),
                    prms => $self->_round($prms),
                    med => $med,
                    min => $min,
                    max => $max,
                    adev => $self->_round($adev),
                    rms => $self->_round($rms),
                );
            }
        }
        if ($self->maximum_depth_position) {
            my ($max_quantity,$max_value) = rle($pdl >= $max);
            my ($max_padded_quantity,$max_padded_value) = $max_quantity->where($max_value,$max_quantity!=0);
            my $max_padded_quantity_cumsum = $max_padded_quantity->cumusumover;
            my ($max_depth_idx,$other_depth_idx) = which_both($max_padded_value);
            my $max_positions = $max_padded_quantity_cumsum->index($max_depth_idx);
            my @positions;
            for (my $i = 0; $i < $max_positions->getdim(0); $i++) {
                my $pdl_position = $max_positions($i)->sclr;
                # I believe the start position is in zero-based coordinates
                my $roi_position = $start + $pdl_position;
                push @positions, $roi_position;
            }
            if ($stats_writer) {
                $data{max_positions} = join(',',@positions);
            }
        }
        if ($stats_writer) {
            if ($data{name}) {
                $stats_writer->write_one(\%data);
            }
        }
    }
    return 1;
}

sub resolve_stats_headers {
    my $self = shift;
    my @headers = @DEFAULT_STATS_HEADERS;
    if ($self->maximum_depth_position) {
        push @headers, 'max_positions';
    }
    return \@headers;
}

1;
