package Genome::Model::Tools::Music::Bmr::CalcCovgHelper;

use warnings;
use strict;
use IO::File;

our $VERSION = $Genome::Model::Tools::Music::VERSION;

class Genome::Model::Tools::Music::Bmr::CalcCovgHelper {
  is => 'Genome::Model::Tools::Music::Bmr::Base',
  has_input => [
    roi_file => { is => 'Text', doc => "Tab delimited list of ROIs [chr start stop gene_name] (See Description)" },
    reference_sequence => { is => 'Text', doc => "Path to reference sequence in FASTA format" },
    normal_tumor_bam_pair => { is => 'Text', doc => "Tab delimited line with sample name, path to normal bam file, and path to tumor bam file (See Description)" },
    output_file => { is => 'Text', doc => "Output file path.  Specify either output-file or output-directory.", is_optional => 1},
    output_dir => { is => 'Text', doc => "Output directory path.  Specify either output-file or output-directory", is_optional => 1},
    normal_min_depth => { is => 'Integer', doc => "The minimum read depth to consider a Normal BAM base as covered", is_optional => 1,  default => 6},
    tumor_min_depth => { is => 'Integer', doc => "The minimum read depth to consider a Tumor BAM base as covered", is_optional => 1, default => 8},
    min_mapq => { is => 'Integer', doc => "The minimum mapping quality of reads to consider towards read depth counts", is_optional => 1, default => 20},
  ],
  has_calculated_optional => [
    sample_name => {
        calculate_from => ['normal_tumor_bam_pair'],
        calculate => q {my @bams = split /\t/, $normal_tumor_bam_pair; return $bams[0];},
    },
    final_output_file => {
        is_output => 1,
        calculate_from => ['output_file','output_dir','sample_name'],
        calculate => q {if ($output_file) {return $output_file;} elsif ($output_dir){return $output_dir."/".$sample_name;} else {die "Either output-file or output-dir must be specified."}},
    },
    normal_bam => {
        calculate_from => ['normal_tumor_bam_pair'],
        calculate => q {my @bams = split /\t/, $normal_tumor_bam_pair; return $bams[1];},
    },
    tumor_bam => {
        calculate_from => ['normal_tumor_bam_pair'],
        calculate => q {my @bams = split /\t/, $normal_tumor_bam_pair; return $bams[2];},
    },
  ],
  doc => "Uses calcRoiCovg.c to count covered bases per-gene for a tumor-normal pair of BAMs."
};

sub help_synopsis {
  return <<HELP
General usage:

 ... music bmr calc-covg-helper \\
    --normal-tumor-bam-pair "sample-name path/to/normal_bam path/to/tumor_bam" \\
    --reference-sequence input_dir/all_sequences.fa \\
    --output-file output_file \\
    --roi-file input_dir/all_coding_exons.tsv

HELP
}

sub help_detail {
  return <<HELP;
This script counts bases with sufficient coverage in the ROIs of each gene in the given pair of
tumor-normal BAM files and categorizes them into - AT, CG (non-CpG), and CpG counts. It also adds
up these base-counts across all ROIs of each gene in the sample, but covered bases that lie
within overlapping ROIs are not counted more than once towards these total counts.

=back
HELP
}

sub _additional_help_sections {
  return (
    "ARGUMENTS",
<<EOS

=over 4

=item --roi-file

=over 8

=item The regions of interest (ROIs) of each gene are typically regions targeted for sequencing or are
  merged exon loci (from multiple transcripts) of genes with 2-bp flanks (splice junctions). ROIs
  from the same chromosome must be listed adjacent to each other in this file. This allows the
  underlying C-based code to run much more efficiently and avoid re-counting bases seen in
  overlapping ROIs (for overall covered base counts). For per-gene base counts, an overlapping
  base will be counted each time it appears in an ROI of the same gene. To avoid this, be sure to
  merge together overlapping ROIs of the same gene. BEDtools' mergeBed can help if used per gene.

=back

=item --reference-sequence

=over 8

=item The reference sequence in FASTA format. If a reference sequence index is not found next to this
  file (a .fai file), it will be created.

=back

=item --normal-tumor-bam-pair

=over 8

=item "sample-name path/to/normal_bam path/to/tumor_bam"

=back

=item --output-file

=over 8

=item Specify an output file where the per-ROI covered base counts will be written

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
  my $tumor_bam = $self->tumor_bam;
  my $normal_bam = $self->normal_bam;
  my $output_file = $self->final_output_file;
  my $normal_min_depth = $self->normal_min_depth;
  my $tumor_min_depth = $self->tumor_min_depth;
  my $min_mapq = $self->min_mapq;

  # Check on all the input data before starting work
  print STDERR "ROI file not found or is empty: $roi_file\n" unless( -s $roi_file );
  print STDERR "Reference sequence file not found: $ref_seq\n" unless( -e $ref_seq );
  print STDERR "Normal BAM file not found or is empty: $normal_bam\n" unless( -s $normal_bam );
  print STDERR "Tumor BAM file not found or is empty: $tumor_bam\n" unless( -s $tumor_bam );
  return undef unless( -s $roi_file && -e $ref_seq && -s $normal_bam && -s $tumor_bam );

  # Check whether the annotated regions of interest are clumped together by chromosome
  my $roiFh = IO::File->new( $roi_file ) or die "ROI file could not be opened. $!\n";
  my @chroms = ( "" );
  while( my $line = $roiFh->getline ) # Emulate Unix's uniq command on the chromosome column
  {
    my ( $chrom ) = ( $line =~ m/^(\S+)/ );
    push( @chroms, $chrom ) if( $chrom ne $chroms[-1] );
  }
  $roiFh->close;
  my %chroms = map { $_ => 1 } @chroms; # Get the actual number of unique chromosomes
  if( scalar( @chroms ) != scalar( keys %chroms ))
  {
    print STDERR "ROIs from the same chromosome must be listed adjacent to each other in file. ";
    print STDERR "If in UNIX, try:\nsort -k 1,1 $roi_file\n";
    return undef;
  }

  # If the reference sequence FASTA file hasn't been indexed, do it
  my $ref_seq_idx = "$ref_seq.fai";
  system( "samtools faidx $ref_seq" ) unless( -e $ref_seq_idx );

  $normal_bam = '' unless( defined $normal_bam );
  $tumor_bam = '' unless( defined $tumor_bam );
  print STDERR "Normal BAM not found: \"$normal_bam\"\n" unless( -e $normal_bam );
  print STDERR "Tumor BAM not found: \"$tumor_bam\"\n" unless( -e $tumor_bam );
  next unless( -e $normal_bam && -e $tumor_bam );

  # Construct the command that calculates coverage per ROI
  my $calcRoiCovg_cmd = "calcRoiCovg $normal_bam $tumor_bam $roi_file $ref_seq $output_file $normal_min_depth $tumor_min_depth $min_mapq";

  # If the calcRoiCovg output was already generated, then don't rerun it
  if( -s $output_file )
  {
    print "Output file $output_file found. Skipping re-calculation.\n";
  }
  # Run the calcRoiCovg command on this tumor-normal pair. This could take a while
  elsif( system( "$calcRoiCovg_cmd" ) != 0 )
  {
    print STDERR "Failed to execute: $calcRoiCovg_cmd\n";
    return;
  }
  else
  {
    print "$output_file generated and stored.\n";
    return 1;
  }

}

1;
