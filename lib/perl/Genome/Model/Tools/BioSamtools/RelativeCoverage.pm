package Genome::Model::Tools::BioSamtools::RelativeCoverage;

use strict;
use warnings;

use Genome;

class Genome::Model::Tools::BioSamtools::RelativeCoverage {
    is => ['Genome::Model::Tools::BioSamtools'],
    has_input => [
        bam_file => {
            doc => 'A BAM file, sorted and indexed, containing alignment/read data',
        },
        bed_file => {
            doc => 'The BED format file (tab delimited: chr,start,end,name) file containing annotation or regions of interest.',
        },
        stats_file => {
            doc => 'When run in parallel, do not define.  From the command line this file will contain the output metrics for each region.',
            is_optional => 1,
        },
        bias_file => {
            doc => 'When run in parallel, do not define.  From the command line this file will contain the output metrics for each region.',
            is_optional => 1,
        },
    ],
};

sub help_detail {
'
These commands are setup to run perl v5.10.0 scripts that use Bio-Samtools and require bioperl v1.6.0.  They all require 64-bit architecture.

Output file format(stats_file):
[1] Region Name (column 4 of BED file)
[2] Percent of Reference Bases Covered
[3] Total Number of Reference Bases
[4] Total Number of Covered Bases
[5] Number of Missing Bases
[6] Average Coverage Depth
[7] Standard Deviation Average Coverage Depth
[8] Median Coverage Depth
[9] Number of Gaps
[10] Average Gap Length
[11] Standard Deviation Average Gap Length
[12] Median Gap Length
[13] Min. Depth Filter
[14] Discarded Bases (Min. Depth Filter)
[15] Percent Discarded Bases (Min. Depth Filter)
';
}


sub execute {
    my $self = shift;
    unless (Genome::Config->arch_os =~ /64/) {
        die('Failed to run on 64-bit architecture');
    }
    my $regions = Genome::Model::Tools::RefCov::ROI::Bed->create(
        file => $self->bed_file,
    );
    unless ($regions) {
        die('Failed to load region file '. $self->bed_file .'.  Accepted formats are: bed');
    }

    open( my $stats_fh, '>'. $self->stats_file ) || die 'Failed to open stats file for writing '. $self->stats_file;

    my $refcov_bam = Genome::Model::Tools::RefCov::Bam->create(
        bam_file => $self->bam_file
    );
    # create low level bam object
    my $bam  = $refcov_bam->bio_db_bam;

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
    my $index = $refcov_bam->bio_db_index;

    my %size_to_relative_depth;
    while (my $region = $regions->next_region) {
        my $id = $region->{name};
        my $target = $region->{chrom};
        my $start = $region->{start};
        my $stop = $region->{end};
        my $length = ($stop - $start) + 1;

        # Here we get the $tid from the $gene_name or $seq_id
        my $tid = $target_name_index{$target};
        unless (defined $tid) { die('Failed to get tid for target '. $target); }

        #low-level API uses zero based coordinates
        #Not sure if this method needs $stop - 1  instead of $stop
        #Needs further testing but thought there was a reason for using just $stop
        my $coverage = $index->coverage( $bam, $tid, $start - 1, $stop );
        unless (scalar( @{ $coverage } ) == $length) {
            die('The length of the ref '. $length .' does not match the depth span '. scalar( @{ $coverage }));
        }
        my $myCoverageStat = Genome::Model::Tools::RefCov::Stats->create( coverage => $coverage);
        print $stats_fh join ("\t", $id, @{ $myCoverageStat->stats_array_ref() }) . "\n";
        if ($self->bias_file) {
            my $size_bin = undef;
            if (($length >= 100) && ($length <= 2_999)) {
                $size_bin = 'SMALL';
            }
            elsif (($length >= 3_000) && ($length <= 6_999)) {
                $size_bin = 'MEDIUM';
            }
            elsif (($length >= 7_000)) {
                $size_bin = 'LARGE';
            }
            my $relative_coverage = Genome::Model::Tools::RefCov::RelativeCoverage->create( coverage => $coverage );
            if (defined($size_bin)) {
                my $hash_ref = $relative_coverage->relative_coverage;
                for my $relative_position (keys %{$hash_ref}) {
                    $size_to_relative_depth{$size_bin}{$relative_position} += $hash_ref->{$relative_position};
                }
            }
        }
    }
    $stats_fh->close;

    if ($self->bias_file) {
        foreach my $size (keys %size_to_relative_depth){
            my $size_bias_file = $self->bias_file .'_'. $size;
            open(my $bias_fh, '>'. $size_bias_file) || die ('Failed to open '. $size .' bias file for writing: '. $size_bias_file);
            my %relative_depth = %{$size_to_relative_depth{$size}};
            my @positions = grep { $_ <= 1 } sort {$a <=> $b} keys %relative_depth;
            my @depth;
            for my $position (@positions) {
                print $bias_fh $position ."\t". $relative_depth{$position} ."\n";
                push @depth, $relative_depth{$position};
            }
            $bias_fh->close;
        }
    }
    return 1;
}

1;
