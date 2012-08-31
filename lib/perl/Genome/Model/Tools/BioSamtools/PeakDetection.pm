package Genome::Model::Tools::BioSamtools::PeakDetection;

use strict;
use warnings;

use Genome;

my $DEFAULT_NORM_FACTORS = '1,2,3,4,5';
my $DEFAULT_OFFSET = 10_000_000;

class Genome::Model::Tools::BioSamtools::PeakDetection {
    is => 'Genome::Model::Tools::BioSamtools',
    has_input => [
        bam_file => {
            is => 'Text',
            doc => 'A path to a BAM format file of aligned capture reads',
        },
        offset => {
            is => 'Integer',
            doc => 'The size(bp) of the offset for the sliding window.',
            default_value => $DEFAULT_OFFSET,
            is_optional => 1,
        },
        normalization_factors => {
            is => 'Text',
            doc => 'A comma separated list of normalization factors to evaluate coverage',
            default_value => $DEFAULT_NORM_FACTORS,
            is_optional => 1,
        },
        output_directory => {
            is => 'Text',
            doc => 'The output directory to generate peak files',
        },
    ],
};

sub execute {
    my $self = shift;
    my $output_directory = $self->output_directory;
    unless (-d $output_directory) {
        unless (Genome::Sys->create_directory($output_directory)) {
            die('Failed to create output_directory: '. $output_directory);
        }
    }
    my $output_file = $output_directory .'/peak_detection.out';
    my $output_fh = Genome::Sys->open_file_for_writing($output_file);
    my @factors = split(',',$self->normalization_factors);
    
    # create low level bam object
    my $refcov_bam  = Genome::Model::Tools::RefCov::Bam->new(bam_file => $self->bam_file );
    unless ($refcov_bam) {
        die('Failed to load bam file '. $self->bam_file);
    }
    my $bam  = $refcov_bam->bio_db_bam;
    my $index = $refcov_bam->bio_db_index;
    my $header = $bam->header();

    # Number of reference sequences
    my $targets = $header->n_targets();

    my $offset = $self->offset;
    my $genome_coverage_stats = Statistics::Descriptive::Sparse->new();
    for (my $tid = 0; $tid < $targets; $tid++) {
        my $chr = $header->target_name->[$tid];
        my $chr_length = $header->target_len->[$tid];
        my %min_depth_previous_clusters;
        for (my $start = 0; $start <= $chr_length; $start += $offset) {
            my $end = $start + $offset;
            if ($end > $chr_length) {
                $end = $chr_length;
            }
            my $coverage = $index->coverage( $bam, $tid, $start, $end);
            my @coverage = @{$coverage};
            my @sorted_coverage = sort {$a <=> $b} (@coverage);
            my $depth = shift(@sorted_coverage);
            while (defined($depth) && $depth == 0) {
                $depth = shift(@sorted_coverage);
            }
            if ($depth) {
                unshift(@sorted_coverage,$depth);
            }
            $genome_coverage_stats->add_data(@sorted_coverage);
        }
    }
    my $mgc = $genome_coverage_stats->mean;
    print $output_fh 'Reference Bases Covered:  '. $genome_coverage_stats->count ."\n";
    print $output_fh 'Mean Reference Coverage:  '.  $mgc ."\n";
    print $output_fh 'Standard Deviation:  '. $genome_coverage_stats->standard_deviation ."\n";
    print $output_fh 'Variance:  '. $genome_coverage_stats->variance ."\n";
    print $output_fh 'Minimum:  '. $genome_coverage_stats->min ."\n";
    print $output_fh 'Maximum:  '. $genome_coverage_stats->max ."\n";
    print $output_fh 'Normalization Factors:  '. join(',',@factors) ."\n";
    my %min_depths;
    for my $factor (sort { $a <=> $b } @factors) {
        my $min_depth = int($factor * $mgc);
        $min_depths{$min_depth} = 1;
    }
    my @min_depths = sort { $a <=> $b } keys %min_depths;
    print $output_fh 'Minimum Depth Filters:  '. join(',', @min_depths);

    my %min_depth_clusters_fhs;
    for my $min_depth (@min_depths) {
        my $min_depth_clusters_file = $output_directory .'/min_depth_'. $min_depth .'_coverage_clusters.bed';
        open (my $min_depth_clusters_fh, '>'. $min_depth_clusters_file ) ||
            die('Failed to open output file:  '. $min_depth_clusters_file);
        $min_depth_clusters_fhs{$min_depth} = $min_depth_clusters_fh;
    }
    for (my $tid = 0; $tid < $targets; $tid++) {
        my $chr = $header->target_name->[$tid];
        my $chr_length = $header->target_len->[$tid];
        my %min_depth_previous_clusters;
        for (my $start = 0; $start <= $chr_length; $start += $offset) {
            my $end = $start + $offset;
            if ($end > $chr_length) {
                $end = $chr_length;
            }
            
            # low-level API uses zero based coordinates
            # all regions should be zero based, but for some reason the correct length is never returned
            # the API must be expecting BED like inputs where the start is zero based and the end is 1-based
            # you can see in the docs for the low-level Bio::DB::BAM::Alignment class that start 'pos' is 0-based,but calend really returns 1-based
            my $coverage = $index->coverage( $bam, $tid, $start, $end);
            for my $min_depth (@min_depths) {
                my @clusters = cluster_coverage($coverage,$min_depth,$start);
                my $min_depth_clusters_fh = $min_depth_clusters_fhs{$min_depth};
                if (@clusters) {
                    if ($min_depth_previous_clusters{$min_depth}) {
                        my $last_cluster = $min_depth_previous_clusters{$min_depth}->[-1];
                        my $first_cluster = $clusters[0];
                        #print $output_fh Data::Dumper::Dumper($last_cluster, $first_cluster);
                        if (($first_cluster->[0] == $last_cluster->[1])) {
                            #Blunt end clusters: merge
                            $last_cluster->[1] = $first_cluster->[1];
                            shift(@clusters);
                        }
                        for my $cluster (@{$min_depth_previous_clusters{$min_depth}}) {
                            print $min_depth_clusters_fh $chr ."\t". $cluster->[0] ."\t". $cluster->[1] ."\n";
                        }
                    }
                    $min_depth_previous_clusters{$min_depth} = \@clusters;
                }
            }
        }
        # One last pass to print remaining clusters to file
        for my $min_depth (keys %min_depth_previous_clusters) {
            my $min_depth_clusters_fh = $min_depth_clusters_fhs{$min_depth};
            for my $cluster (@{$min_depth_previous_clusters{$min_depth}}) {
                print $min_depth_clusters_fh $chr ."\t". $cluster->[0] ."\t". $cluster->[1] ."\n";
            }
            $min_depth_previous_clusters{$min_depth} = [];
        }
    }
    for my $min_depth (keys %min_depth_clusters_fhs) {
        my $fh = $min_depth_clusters_fhs{$min_depth};
        $fh->close;
    }
    $output_fh->close;
    return 1;
}

sub cluster_coverage {
    my ($coverage, $min_depth, $offset) = @_;
    my @clusters;
    # Populate the cluster data structure based on the coverage depth array.

    map {$_ = undef} my ($is_cluster, $start, $stop);
    my $last = scalar(@{ $coverage });
    for (my $i = 0; $i < $last; $i++) {
        my $depth = $coverage->[$i];
        if ($depth >= $min_depth) {
            if ($is_cluster) {
                # continue extending existing cluster
                if ($i == ($last - 1)) {
                    # end of coverage array, so end cluster and finish
                    $stop = ($i + 1) + $offset;
                } else {
                    #reading through coverage
                    next;
                }
            } else {
                # start new cluster
                unless ($offset) {
                    $start = 1;
                }
                $start = $i + $offset;
                $is_cluster = 1;
            }
        } else {
            if ($is_cluster) {
                # end of cluster, start of gap
                $stop = $i + $offset;
                $is_cluster = undef;
            } else {
                # reading through a gap
                next;
            }
        }
        # Update collection hash:
        if (defined($start) && defined($stop)) {
            push @clusters, [$start,$stop];
            $start = undef;
            $stop = undef;
        } elsif ($i == ($last - 1)) {
            # one base long cluster on end of coverage array
            push @clusters, [$start,$start];
            $start = undef;
            $stop = undef;
        }
    }
    return @clusters;
}


1;
