# Author: Cyriac Kandoth
# Date: 11.07.2010
# Description: Creates a wiggle track format file for a tumor-normal pair of BAM files (uses bam2wig)

package Genome::Model::Tools::Capture::WiggleFromBams;
use strict;
use warnings;
use IO::File;
use Benchmark;

class Genome::Model::Tools::Capture::WiggleFromBams {
    is => 'Command',
    has => [
        normal_bam => { is => 'File', doc => "Normal BAM", is_optional => 0 },
        tumor_bam => { is => 'File', doc => "Tumor BAM", is_optional => 0 },
        regions_file => { is => 'File', doc => "Targeted regions: 3-column, tab-delimited, 1-based loci", is_optional => 0 },
        output_file => { is => 'File', doc => "Wiggle file to output (.wig extension recommended)", is_optional => 0 },
        min_base_qual => { is => 'Integer', doc => "Minimum base quality of reads to count", is_optional => 1, default => 20 },
        min_depth_normal => { is => 'Integer', doc => "Minimum Q20 depth for Normal", is_optional => 1, default => 6 },
        min_depth_tumor => { is => 'Integer', doc => "Minimum Q20 depth for Tumor", is_optional => 1, default => 8 },
    ],
};

sub help_brief {
    return "Generates a wiggle file given Tumor/Normal BAMs and the targeted regions";
}

sub help_detail {
    return <<EOS
Builds a wiggle track format file using fixedStep=1, for a given list of targeted regions,
and two BAM files (tumor and normal). Two intermediate normal and tumor wiggles files are created
with the extensions wig_normal and wig_tumor. They can be removed if necessary.
EOS
}

sub execute {
    my $self = shift;
    my $t0 = Benchmark->new;

    my $normal_bam = $self->normal_bam;
    my $tumor_bam = $self->tumor_bam;
    my $regions_file = $self->regions_file;
    my $output_wig = $self->output_file;
    my $min_base_qual = $self->min_base_qual;
    my $min_depth_normal = $self->min_depth_normal;
    my $min_depth_tumor = $self->min_depth_tumor;

    print "Checking BAM files...\n";
    $self->error( "Normal BAM not found ($normal_bam)." ) unless( -e $normal_bam );
    $self->error( "Tumor BAM not found ($tumor_bam)." ) unless( -e $tumor_bam );

    print "Checking BAM index files...\n";
    unless( -e "$normal_bam.bai" ) {
        print "Normal BAM index file not found. Using \"samtools index\" to create it...\n";
        my $stdoe = `samtools index $normal_bam 2>&1`;
        $self->error( "samtools failed! Be sure to use 64-bit architecture." ) if( $stdoe =~ m/syntax error: \S+ unexpected/i );
    }
    unless( -e "$tumor_bam.bai" ) {
        print "Tumor BAM index file not found. Using \"samtools index\" to create it...\n";
        my $stdoe = `samtools index $tumor_bam 2>&1`;
        $self->error( "samtools failed! Be sure to use 64-bit architecture." ) if( $stdoe =~ m/syntax error: \S+ unexpected/i );
    }

    print "Checking targeted regions file...\n";
    my $inRgnFh = IO::File->new( $regions_file );
    my $line = $inRgnFh->getline; # Potentially a header. Don't allow that
    $self->error( "Cannot recognize format of Line 1 in the targeted regions file." ) if( $line !~ m/^\w+\t\d+\t\d+/ );
    $line = $inRgnFh->getline while( $line !~ m/^\w+\t\d+\t\d+/ ); # Find a line that's actually data
    $self->error( "Please remove 'chr' prefixes for chromosome names in the targeted regions file." ) if( $line =~ m/^chr/i );
    $inRgnFh->close;

    print "Using bam2wig to generate wiggle files for tumor and normal BAMs...\n";
    $self->error( "bam2wig failed! Be sure to use 64-bit architecture." ) if( `bam2wig 2>&1` =~ m/syntax error: \S+ unexpected/i );
    my ( $normal_wig, $tumor_wig ) = ( "$output_wig\_normal", "$output_wig\_tumor" );
    `bam2wig -q $min_base_qual -d $min_depth_normal -o $normal_wig -l $regions_file $normal_bam`;
    `bam2wig -q $min_base_qual -d $min_depth_tumor -o $tumor_wig -l $regions_file $tumor_bam`;

    print "Merging the two wiggle files into one consolidated wiggle file...\n";
    my $outWigFh = IO::File->new( $output_wig, ">" );
    my $inNWigFh = IO::File->new( $normal_wig );
    my $inTWigFh = IO::File->new( $tumor_wig );
    while( my $n_line = $inNWigFh->getline ) {
        my $t_line = $inTWigFh->getline;
        if( $n_line !~ m/^1$/ ) {
            #prints fixedStep lines as well as zeroes
            $outWigFh->print( $n_line );
        }
        else {
            $outWigFh->print( $t_line );
        }
    }
    $inTWigFh->close;
    $inNWigFh->close;
    $outWigFh->close;

    my $t1 = Benchmark->new;
    print "Done! Wall Time: ", timestr( timediff( $t1, $t0 )), "\n";
    return 1;
}

1;
