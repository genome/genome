package Genome::Model::Tools::Music::Bmr::CalcWigCovg;

use warnings;
use strict;
use Genome;
use IO::File;

our $VERSION = $Genome::Model::Tools::Music::VERSION;

class Genome::Model::Tools::Music::Bmr::CalcWigCovg {
    is => 'Genome::Model::Tools::Music::Bmr::Base',
    has_input => [
        roi_file => { is => 'Text', doc => "Tab-delimited list of ROIs [chr start stop gene_name] (See Description)" },
        reference_sequence => { is => 'Text', doc => "Path to reference sequence in FASTA format" },
        wig_list => { is => 'Text', doc => "Tab-delimited list of WIG files [sample_name wig_file] (See Description)" },
        output_dir => { is => 'Text', doc => "Directory where output files and subdirectories will be written", is_output => 1 },
    ],
    doc => "Count covered bases per-gene for each given wiggle track format file."
};

sub help_synopsis {
    return <<HELP
General usage:

 ... music bmr calc-wig-covg \\
    --wig-list input_dir/wig_list \\
    --output-dir output_dir/ \\
    --reference-sequence input_dir/all_sequences.fa \\
    --roi-file input_dir/all_coding_exons.tsv
HELP
}

sub help_detail {
    return <<HELP;
This script counts bases with sufficient coverage in the ROIs of each gene from given wiggle track
format files, and categorizes them into - AT, CG (non-CpG), and CpG counts. It also adds
up these base-counts across all ROIs of each gene for each sample, but covered bases that lie
within overlapping ROIs are not counted more than once towards these total counts.
HELP
}

sub _additional_help_sections {
    return (
    "ARGUMENTS",
<<EOS

=over 4

=item --roi-file

=over 8

=item The regions of interest (ROIs) of each gene are typically regions targeted for sequencing or
  are merged exon loci (from multiple transcripts) of genes with 2-bp flanks (splice junctions).
  For per-gene base counts, an overlapping base will be counted each time it appears in an ROI of
  the same gene. To avoid this, be sure to merge together overlapping ROIs of the same gene.
  BEDtools' mergeBed can help if used per gene.

=back

=item --reference-sequence

=over 8

=item The reference sequence in FASTA format. If a reference sequence index is not found next to
  this file (a .fai file), it will be created.

=back

=item --wig-list

=over 8

=item Provide a file containing sample names and the wiggle track format file locations for each.
  Use the tab-delimited format [sample_name wig_file] per line. Additional columns like clinical
  data are allowed, but ignored. The sample_name must be the same as the tumor sample names used
  in the MAF file (16th column, with the header Tumor_Sample_Barcode).

=back

=item --output-dir

=over 8

=item Specify an output directory where the following will be created/written:
  roi_covgs: Subdirectory containing per-ROI covered base counts for each sample.
  gene_covgs: Subdirectory containing per-gene covered base counts for each sample.
  total_covgs: File containing the overall non-overlapping coverages per sample.

=back

=back

EOS
    );
}

sub _doc_authors {
    return " Cyriac Kandoth, Ph.D.";
}

sub _doc_see_also {
    return <<EOS
B<genome-music-bmr>(1),
B<genome-music>(1),
B<genome>(1)
EOS
}

sub execute {
    my $self = shift;
    my $roi_file = $self->roi_file;
    my $ref_seq = $self->reference_sequence;
    my $wig_list = $self->wig_list;
    my $output_dir = $self->output_dir;

    # Check on all the input data before starting work
    print STDERR "ROI file not found or is empty: $roi_file\n" unless( -s $roi_file );
    print STDERR "Reference sequence file not found: $ref_seq\n" unless( -e $ref_seq );
    print STDERR "List of WIGs not found or is empty: $wig_list\n" unless( -s $wig_list );
    print STDERR "Output directory not found: $output_dir\n" unless( -e $output_dir );
    return undef unless( -s $roi_file && -e $ref_seq && -s $wig_list && -e $output_dir );

    # Outputs of this script will be written to these locations in the output directory
    $output_dir =~ s/(\/)+$//; # Remove trailing forward slashes if any
    my $roi_covg_dir = "$output_dir/roi_covgs"; # Stores output from calcRoiCovg per sample
    my $gene_covg_dir = "$output_dir/gene_covgs"; # Stores per-gene coverages per sample
    my $tot_covg_file = "$output_dir/total_covgs"; # Stores total coverages per sample

    # If the reference sequence FASTA file hasn't been indexed, do it
    my $ref_seq_idx = "$ref_seq.fai";
    unless( -e $ref_seq_idx ) {
        print "Reference fasta index not found. Creating one at: $ref_seq.fai\n";
        system( "samtools faidx $ref_seq" ) or die "Failed to run samtools! $!\n";
    }

    # Create a temporary 0-based ROI BED-file that we can use with joinx, and also measure gene lengths
    my %geneLen = ();
    my $roi_bed = Genome::Sys->create_temp_file_path();
    my $roiBedFh = IO::File->new( $roi_bed, ">" ) or die "Temporary ROI BED file could not be created. $!\n";
    my $roiFh = IO::File->new( $roi_file ) or die "ROI file could not be opened. $!\n";
    while( my $line = $roiFh->getline ) {
        chomp( $line );
        my ( $chr, $start, $stop, $gene ) = split( /\t/, $line );
        --$start;
        unless( $start >= 0 && $start < $stop ) {
            print STDERR "Invalid ROI: $line\nPlease use 1-based loci and ensure that start <= stop\n";
            return undef;
        }
        $geneLen{$gene} += ( $stop - $start );
        $roiBedFh->print( "$chr\t$start\t$stop\t$gene\n" );
    }
    $roiFh->close;
    $roiBedFh->close;

    # Also create a merged BED file where overlapping ROIs are joined together into contiguous regions
    # ::TODO:: Use joinx instead of mergeBed, because we'd rather add an in-house dependency
    my $merged_roi_bed = Genome::Sys->create_temp_file_path();
    system( "mergeBed -i $roi_bed | joinx1.7 sort -s - -o $merged_roi_bed" );# or die "Failed to run mergeBed or joinx!\n$roi_bed\n$merged_roi_bed\n $!\n";

    # Create the output directories unless they already exist
    mkdir $roi_covg_dir unless( -e $roi_covg_dir );
    mkdir $gene_covg_dir unless( -e $gene_covg_dir );

    # This is a file that will report the overall non-overlapping coverages per WIG
    my $totCovgFh = IO::File->new( $tot_covg_file, ">" );
    $totCovgFh->print( "#Sample\tCovered_Bases\tAT_Bases_Covered\tCG_Bases_Covered\tCpG_Bases_Covered\n" );

    # Parse through each pair of WIG files provided and run calcRoiCovg as necessary
    my $wigFh = IO::File->new( $wig_list );
    while( my $line = $wigFh->getline ) {
        next if( $line =~ m/^#/ );
        chomp( $line );
        my ( $sample, $wig_file ) = split( /\t/, $line );
        $wig_file = '' unless( defined $wig_file );
        print STDERR "Wiggle track format file for $sample not found: \"$wig_file\"\n" unless( -e $wig_file );
        next unless( -e $wig_file );

        # Use joinx to parse the WIG file and return per-ROI coverages of AT, CG (non-CpG), and CpG
        system( "joinx1.7 wig2bed -Zc $wig_file | joinx1.7 sort -s | joinx1.7 intersect -F \"I A3\" $roi_bed - | joinx1.7 ref-stats - $ref_seq | cut -f 1-7 > $roi_covg_dir/$sample.covg" );# or die "Failed to run joinx to calculate per-gene coverages in $sample! $!\n";

        # Read the joinx formatted coverage file and count covered bases per gene
        my %geneCovg = ();
        my $roiCovgFh = IO::File->new( "$roi_covg_dir/$sample.covg" );
        while( my $line = $roiCovgFh->getline ) {
            chomp( $line );
            if( $line !~ m/^#/ ) {
                my ( undef, undef, undef, $gene, $at_covd, $cg_covd, $cpg_covd ) = split( /\t/, $line );
                $geneCovg{$gene}{covd} += ( $at_covd + $cg_covd + $cpg_covd );
                $geneCovg{$gene}{at} += $at_covd;
                $geneCovg{$gene}{cg} += $cg_covd;
                $geneCovg{$gene}{cpg} += $cpg_covd;
            }
        }
        $roiCovgFh->close;

        # Write the per-gene coverages to a file named after this sample_name
        my $geneCovgFh = IO::File->new( "$gene_covg_dir/$sample.covg", ">" );
        $geneCovgFh->print( "#Gene\tLength\tCovered\tAT_covd\tCG_covd\tCpG_covd\n" );
        foreach my $gene ( sort keys %geneLen ) {
            if( defined $geneCovg{$gene} ) {
                $geneCovgFh->print( join( "\t", $gene, $geneLen{$gene}, $geneCovg{$gene}{covd},
                                    $geneCovg{$gene}{at}, $geneCovg{$gene}{cg}, $geneCovg{$gene}{cpg} ), "\n" );
            }
            else {
                $geneCovgFh->print( "$gene\t" . $geneLen{$gene} . "\t0\t0\t0\t0\n" );
            }
        }
        $geneCovgFh->close;

        # Measure coverage stats on the merged ROI file, so that bps across the genome are not counted twice
        my $merged_roi_bed_covg = Genome::Sys->create_temp_file_path();
        system( "joinx1.7 wig2bed -Zc $wig_file | joinx1.7 sort -s | joinx1.7 intersect $merged_roi_bed - | joinx1.7 ref-stats - $ref_seq | cut -f 1-6 > $merged_roi_bed_covg" );# or die "Failed to run joinx to calculate overall coverages in $sample! $!\n";

        # Read the joinx formatted coverage file and sum up the coverage stats per region
        my ( $tot_covd, $tot_at_covd, $tot_cg_covg, $tot_cpg_covd );
        my $totRoiCovgFh = IO::File->new( $merged_roi_bed_covg );
        while( my $line = $totRoiCovgFh->getline ) {
            chomp( $line );
            if( $line !~ m/^#/ ) {
                my ( $chr, $start, $stop, $at_covd, $cg_covd, $cpg_covd ) = split( /\t/, $line );
                $tot_covd += ( $at_covd + $cg_covd + $cpg_covd );
                $tot_at_covd += $at_covd;
                $tot_cg_covg += $cg_covd;
                $tot_cpg_covd += $cpg_covd;
            }
        }
        $totRoiCovgFh->close;
        $totCovgFh->print( "$sample\t$tot_covd\t$tot_at_covd\t$tot_cg_covg\t$tot_cpg_covd\n" );
    }
    $wigFh->close;
    $totCovgFh->close;

    return 1;
}

1;
