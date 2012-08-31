package Genome::Model::Tools::BioSamtools::RefCovCompare;

use strict;
use warnings;

use Genome;

class Genome::Model::Tools::BioSamtools::RefCovCompare {
    is => ['Genome::Model::Tools::BioSamtools'],
    has_input => [
        bam_file_a => {
            doc => 'A BAM file, sorted and indexed, containing alignment/read data',
        },
        bam_file_b => {
            doc => 'A BAM file, sorted and indexed, containing alignment/read data',
        },
        bed_file => {
            doc => 'The BED format file (tab delimited: chr,start,end,name) file containing annotation or regions of interest.',
        },
        output_file => {
            doc => 'The path where an output file will be written with A<->B coverage comparison',
        },
        min_depth_filter => {
            doc => 'The minimum depth at each position to consider coverage.',
            default_value => 1,
            is_optional => 1,
        },
    ],
};

sub help_detail {
'
These commands are setup to run perl v5.10.0 scripts that use Bio-Samtools and require bioperl v1.6.0.  They all require 64-bit architecture.

Output file format:
[1] Region Name (column 4 of BED file)
[2] Region Length
[3] AB - Missing Base Pair
[4] AB - Percent Missing Base Pair
[5] AB - Total Covered Base Pair
[6] AB - Percent Total Covered Base Pair
[7] AB - Unique Covered Base Pair
[8] AB - Perecent Unique Covered Base Pair
[9] A - Total Covered Base Pair
[10] A - Percent Total Covered Base Pair
[11] A - Unique Covered Base Pair
[12] A - Percent Unique Covered Base Pair
[13] B - Total Covered Base Pair
[14] B - Percent Total Covered Base Pair
[15] B - Unique Covered Base Pair
[16] B - Percent Unique Covered Base Pair

'
}


sub execute {
    my $self = shift;

    my $regions = Genome::Model::Tools::RefCov::ROI::Bed->create(
        file => $self->bed_file
    );
    unless ($regions) {
        die('Failed to load region file '. $self->bed_file .'.  Accepted formats are: bed');
    }
    my $min_depth = $self->min_depth_filter;
    unless ($min_depth) {
        $min_depth = 1;
    }
    my $out_fh = IO::File->new($self->output_file,'w');
    unless ($out_fh) {
        die('Failed to open file for output '. $self->output_file);
    }

    # create low level bam object
    my $refcov_bam_a  = Genome::Model::Tools::RefCov::Bam->create(bam_file => $self->bam_file_a );
    my $refcov_bam_b  = Genome::Model::Tools::RefCov::Bam->create(bam_file => $self->bam_file_b );

    my $bam_a = $refcov_bam_a->bio_db_bam;
    my $bam_b = $refcov_bam_b->bio_db_bam;

    my $header_a = $bam_a->header();
    my $header_b = $bam_b->header();

    # Number of reference sequences
    my $targets_a = $header_a->n_targets();
    my $targets_b = $header_b->n_targets();
    unless ($targets_a eq $targets_b) {
        die('BAM file a '. $self->bam_file_a .' had '. $targets_a .' targets and BAM file b '. $self->bam_file_b .' had '. $targets_b .' targets.');
    }

    # The reference sequence names in an array ref with indexed positions
    my $target_names_a = $header_a->target_name();
    my $target_names_b = $header_b->target_name();


    # at the low level API the seq_id/target_name is meaningless
    # cache the target_names in a hash by actual reference sequence name
    # then we can look up the target index on the fly
    my %target_name_index;
    my $i;
    for ( $i = 0; $i < $targets_a; $i++) {
        my $target_name_a = $target_names_a->[$i];
        my $target_name_b = $target_names_b->[$i];
        unless ($target_name_a eq $target_name_b) {
            die('The target names are not the same or out of order.  The '. $i .' target in file A '. $target_name_a .' and in file B '. $target_name_b);
        }
        $target_name_index{$target_name_a} = $i;
    }

    # Make sure our index is not off
    unless ($targets_a == $i) {
        die 'Expected '. $targets_a .' targets but counted '. $i .' indices';
    }
    
    my $index_a = $refcov_bam_a->bio_db_index;
    my $index_b = $refcov_bam_b->bio_db_index;
    
    #Outdated header
    #print "ID\tUNCOVERED\tUNIQUE:$bam_file_a\tSHARED\tUNIQUE:$bam_file_b\n";
    while ( my $region = $regions->next_region) {
        my $id = $region->{name};
        my $target = $region->{chrom};
        $target =~ s/chr//;
        my $start = $region->{start};
        my $stop = $region->{end};
        my $length = ($stop - $start) + 1;
        # Here we get the $tid from the $gene_name or $seq_id
        my $tid = $target_name_index{$target};
        unless (defined $tid) { die('Failed to get tid for target '. $target); }

        #low-level API uses zero based coordinates
        #Not sure if this method needs $stop - 1  instead of $stop
        #Needs further testing but thought there was a reason for using just $stop
        my $coverage_a = $index_a->coverage( $bam_a, $tid, $start - 1, $stop );
        unless (scalar( @{ $coverage_a } ) == $length) {
            die('The length of the locus '. $length .' does not match the coverage length '. scalar( @{ $coverage_a }));
        }
        my $coverage_b = $index_b->coverage( $bam_b, $tid, $start - 1, $stop );
        unless (scalar( @{ $coverage_b } ) == $length) {
            die('The length of the locus '. $length .' does not match the coverage length '. scalar( @{ $coverage_b }));
        }
        my $refcov_stats_a = Genome::Model::Tools::RefCov::Stats->create( coverage => $coverage_a , min_depth => $min_depth);
        my $refcov_stats_b = Genome::Model::Tools::RefCov::Stats->create( coverage => $coverage_b , min_depth => $min_depth);
        #my ($uncovered, $a_unique, $shared, $b_unique) = &compare_coverage_arrays($coverage_a, $coverage_b);
        my $overlap_results = &overlap_results($coverage_a,$coverage_b);
        #unless ($overlap_results->{no_set_coverage} eq $uncovered) {
        #    die;
        #}
        #print $id ."\t". $uncovered ."\t". $a_unique ."\t". $shared ."\t". $b_unique ."\n";
        #print Data::Dumper::Dumper($overlap_results);
        print $out_fh $id ."\t". $overlap_results->{array_length} ."\t".
            $overlap_results->{no_set_coverage} ."\t". $overlap_results->{no_set_coverage_percent} ."\t".
                $overlap_results->{AB_total_coverage} ."\t". $overlap_results->{AB_total_coverage_percent} ."\t".
                    $overlap_results->{AB_unique_coverage} ."\t". $overlap_results->{AB_unique_coverage_percent} ."\t".
                        $overlap_results->{A_total_coverage} ."\t". $overlap_results->{A_total_coverage_percent} ."\t".
                            $overlap_results->{A_unique_coverage} ."\t". $overlap_results->{A_unique_coverage_percent} ."\t".
                                $overlap_results->{B_total_coverage} ."\t". $overlap_results->{B_total_coverage_percent} ."\t".
                                    $overlap_results->{B_unique_coverage} ."\t". $overlap_results->{B_unique_coverage_percent} ."\n";
    }
    $out_fh->close;
    return 1;
}

sub compare_coverage_arrays {
    my $a_ref = shift;
    my $b_ref = shift;
    unless (scalar(@{$a_ref}) == scalar(@{$b_ref})) {
        die ('Passed array refs with uneven size');
    }
    my $a_unique = 0;
    my $shared = 0;
    my $b_unique = 0;
    my $uncovered = 0;
    my $length = scalar(@{$a_ref});
    for (my $i = 0; $i < $length; $i++) {
        if ($a_ref->[$i] && $b_ref->[$i]) {
            $shared++;
        } elsif ( $a_ref->[$i] ) {
            $a_unique++;
        } elsif ( $b_ref->[$i] ) {
            $b_unique++;
        } else {
            $uncovered++;
        }
    }
    return ($uncovered, $a_unique, $shared, $b_unique);
}


sub overlap_results {
    my ($aref_A, $aref_B) = @_;
    unless (scalar(@{$aref_A}) == scalar(@{$aref_B})) {
        die ('Passed array refs with uneven size');
    }
    # Assumes both arrays are of same length.
    my $array_length = scalar( @{ $aref_A } );
    my $end          = $array_length - 1;
    my @AB;
    for (0..$end) {
        my $AB_pos = 0;
        if (${$aref_A}[$_] >= 1) { $AB_pos++ }
        if (${$aref_B}[$_] >= 1) { $AB_pos++ }
        $AB[$_] = $AB_pos;
    }

    # Re-visit arrays A & B, determine overlaps, etc.
    map {$_ = 0} my (
                     $no_set_coverage,
                     $both_set_coverage,
                     $A_unique_coverage,
                     $B_unique_coverage,
                    );
    for (0..$end) {
        if ($AB[$_] == 0) {
            # No coverage.
            $no_set_coverage++;
        }
        elsif ($AB[$_] == 1) {
            # Single set coverage.
            if (${$aref_A}[$_] >= 1) { $A_unique_coverage++ }
            if (${$aref_B}[$_] >= 1) { $B_unique_coverage++ }
        }
        elsif ($AB[$_] == 2) {
            # Both set coverage.
            $both_set_coverage++;
        }
    }
    my $A_total_coverage = $A_unique_coverage + $both_set_coverage;
    my $B_total_coverage = $B_unique_coverage + $both_set_coverage;
    my $AB_total_coverage = $A_unique_coverage + $B_unique_coverage + $both_set_coverage;
    
    # Update the report hash.
    my $overlap_results = {
                           array_length              => $array_length,
                           no_set_coverage           => $no_set_coverage,
                           AB_total_coverage         => $AB_total_coverage,
                           AB_unique_coverage        => $both_set_coverage,
                           A_total_coverage          => $A_total_coverage,
                           A_unique_coverage         => $A_unique_coverage,
                           B_total_coverage          => $B_total_coverage,
                           B_unique_coverage         => $B_unique_coverage,
                           no_set_coverage_percent   => _round( ($no_set_coverage   / $array_length) * 100 ),
                           AB_total_coverage_percent => _round( ($AB_total_coverage / $array_length) * 100 ),
                           AB_unique_coverage_percent=> _round( ($both_set_coverage / $array_length) * 100 ),
                           A_total_coverage_percent  => _round( ($A_total_coverage / $array_length) * 100 ),
                           A_unique_coverage_percent => _round( ($A_unique_coverage / $array_length) * 100 ),
                           B_total_coverage_percent  => _round( ($B_total_coverage / $array_length) * 100 ),
                           B_unique_coverage_percent => _round( ($B_unique_coverage / $array_length) * 100 ),
                          };


    sub _round { return sprintf( "%.2f", shift ) }

    return $overlap_results;
}


1;
