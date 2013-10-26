package Genome::Model::Tools::Music::Bmr::CalcCovg;

use warnings;
use strict;
use IO::File;

our $VERSION = $Genome::Model::Tools::Music::VERSION;

class Genome::Model::Tools::Music::Bmr::CalcCovg {
  is => 'Genome::Model::Tools::Music::Bmr::Base',
  has_input => [
    roi_file => { is => 'Text', doc => "Tab delimited list of ROIs [chr start stop gene_name] (See Description)" },
    reference_sequence => { is => 'Text', doc => "Path to reference sequence in FASTA format" },
    bam_list => { is => 'Text', doc => "Tab delimited list of BAM files [sample_name normal_bam tumor_bam] (See Description)" },
    output_dir => { is => 'Text', doc => "Directory where output files and subdirectories will be written", is_output => 1},
    cmd_list_file => { is => 'Text', doc => "A file to write calcRoiCovg commands to (See Description)", is_optional => 1 },
    cmd_prefix => { is => 'Text', doc => "A command that submits a job to your cluster (See Description)", is_optional => 1 },
    normal_min_depth => { is => 'Integer', doc => "The minimum read depth to consider a Normal BAM base as covered", is_optional => 1},
    tumor_min_depth => { is => 'Integer', doc => "The minimum read depth to consider a Tumor BAM base as covered", is_optional => 1},
    min_mapq => { is => 'Integer', doc => "The minimum mapping quality of reads to consider towards read depth counts", is_optional => 1},
  ],
  has_output => [
    gene_covg_dir => {is_optional => 1, is => 'Text', doc => "Directory where per-sample gene coverage files are located"},
  ],
  doc => "Uses calcRoiCovg.c to count covered bases per-gene for each given tumor-normal pair of BAMs."
};

sub help_synopsis {
  return <<HELP
General usage:

 ... music bmr calc-covg \\
    --bam-list input_dir/bam_list \\
    --output-dir output_dir/ \\
    --reference-sequence input_dir/all_sequences.fa \\
    --roi-file input_dir/all_coding_exons.tsv

To create a list of commands that will allow the processing of each tumor-normal pair in parallel
with an LSF job scheduler:

 ... music bmr calc-covg \\
    --bam-list input_dir/bam_list \\
    --output-dir output_dir/ \\
    --reference-sequence input_dir/all_sequences.fa \\
    --roi-file input_dir/all_coding_exons.tsv \\
    --cmd_list_file parallelizable_commands \\
    --cmd_prefix bsub

In the above case, the commands printed into the output file "parallelizable_commands" can be run
in parallel. After they complete, rerun this script as printed directly below (--cmd_list_file
and --cmd_prefix have been removed) to merge the parallelized calculations:

 ... music bmr calc-covg \\
    --bam-list input_dir/bam_list \\
    --output-dir output_dir/ \\
    --reference-sequence input_dir/all_sequences.fa \\
    --roi-file input_dir/all_coding_exons.tsv
HELP
}

sub help_detail {
  return <<HELP;
This script counts bases with sufficient coverage in the ROIs of each gene in the given pairs of
tumor-normal BAM files and categorizes them into - AT, CG (non-CpG), and CpG counts. It also adds
up these base-counts across all ROIs of each gene for each sample, but covered bases that lie
within overlapping ROIs are not counted more than once towards these total counts.

By default, this script runs a C-based tool named calcRoiCovg for each sample one after another,
taking ~30 mins per sample to generate per-ROI covered base counts. If the results of calcRoiCovg
for a sample already exists in the output subdirectory roi_covgs, re-calculation is skipped. This
allows you to run your own calcRoiCovg jobs in parallel or on multiple machines (Keep reading).

Speed things up by running calcRoiCovg jobs in parallel:
If a compute cluster or multiple machines are available, run this script twice as follows:

=over 4

=item * Define cmd-list-file and cmd-prefix to generate a file with commands that can be submitted to a
cluster or run manually. These jobs will write per-ROI base counts in a subdirectory roi_covgs.

=item * After all the parallelized calcRoiCovg jobs are completed, run this script again to add them up
and generate the final per-gene base counts in a subdirectory gene_covgs. Remember to remove the
cmd-list-file and cmd-prefix arguments or you will just be re-creating a list of commands.

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

=item --bam-list

=over 8

=item Provide a file containing sample names and normal/tumor BAM locations for each. Use the tab-
  delimited format [sample_name normal_bam tumor_bam] per line. Additional columns like clinical
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

=item --cmd-list-file

=over 8

=item Specify a file into which a list of calcRoiCovg jobs will be written to. These can be scheduled
  in parallel, and will write per-ROI covered base-counts into the output subdirectory roi_covgs.
  If cmd-list-file is left unspecified, this script runs calcRoiCovg per sample one after another,
  taking ~30 mins per sample, but it skips samples whose output is already in roi_covgs.

=back

=item --cmd-prefix

=over 8

=item Specify a job submission command that will be prefixed to each command in cmd-list-file. This
  makes batch submission easier. Just run the cmd-list-file file as a shell script to submit jobs.
  cmd-prefix is "bsub" if your cluster uses the LSF job scheduler, or "qsub" in Torque. Add
  arguments as necessary. For example, "bsub -M 4GB" sets a soft memory limit of 4GB.

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
  my $bam_list = $self->bam_list;
  my $output_dir = $self->output_dir;
  my $cmd_list_file = $self->cmd_list_file;
  my $cmd_prefix = $self->cmd_prefix;
  my $normal_min_depth = $self->normal_min_depth;
  my $tumor_min_depth = $self->tumor_min_depth;
  my $min_mapq = $self->min_mapq;

  my $optional_params = "";

  if ($normal_min_depth) {
    $optional_params .= " --normal-min-depth $normal_min_depth";
  }
  if ($tumor_min_depth) {
    $optional_params .= " --tumor-min-depth $tumor_min_depth";
  }
  if ($min_mapq) {
    $optional_params .= " --min-mapq $min_mapq";
  }

  # Check on all the input data before starting work
  print STDERR "ROI file not found or is empty: $roi_file\n" unless( -s $roi_file );
  print STDERR "Reference sequence file not found: $ref_seq\n" unless( -e $ref_seq );
  print STDERR "List of BAMs not found or is empty: $bam_list\n" unless( -s $bam_list );
  print STDERR "Output directory not found: $output_dir\n" unless( -e $output_dir );
  return undef unless( -s $roi_file && -e $ref_seq && -s $bam_list && -e $output_dir );

  # Outputs of this script will be written to these locations in the output directory
  $output_dir =~ s/(\/)+$//; # Remove trailing forward slashes if any
  my $roi_covg_dir = "$output_dir/roi_covgs"; # Stores output from calcRoiCovg per sample
  my $gene_covg_dir = "$output_dir/gene_covgs"; # Stores per-gene coverages per sample
  my $tot_covg_file = "$output_dir/total_covgs"; # Stores total coverages per sample

  $self->gene_covg_dir($gene_covg_dir);

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

  # Create the output directories unless they already exist
  mkdir $roi_covg_dir unless( -e $roi_covg_dir );
  mkdir $gene_covg_dir unless( -e $gene_covg_dir );

  my ( $cmdFh, $totCovgFh );
  if( defined $cmd_list_file )
  {
    $cmdFh = IO::File->new( $cmd_list_file, ">" );
    print "Creating a list of parallelizable jobs at $cmd_list_file.\n";
    print "After successfully running all the jobs in $cmd_list_file,\n",
          "be sure to run this script a second time (without defining the cmd-list-file argument) to merge results in roi_covgs.\n";
  }
  else
  {
    $totCovgFh = IO::File->new( $tot_covg_file, ">" );
    $totCovgFh->print( "#Sample\tCovered_Bases\tAT_Bases_Covered\tCG_Bases_Covered\tCpG_Bases_Covered\n" );
  }

  # Parse through each pair of BAM files provided and run calcRoiCovg as necessary
  my $bamFh = IO::File->new( $bam_list );
  while( my $line = $bamFh->getline )
  {
    next if( $line =~ m/^#/ );
    chomp( $line );
    my ( $sample, $normal_bam, $tumor_bam ) = split( /\t/, $line );
    $normal_bam = '' unless( defined $normal_bam );
    $tumor_bam = '' unless( defined $tumor_bam );
    print STDERR "Normal BAM for $sample not found: \"$normal_bam\"\n" unless( -e $normal_bam );
    print STDERR "Tumor BAM for $sample not found: \"$tumor_bam\"\n" unless( -e $tumor_bam );
    next unless( -e $normal_bam && -e $tumor_bam );

    # Construct the command that calculates coverage per ROI
    my $calcRoiCovg_cmd = "\'gmt music bmr calc-covg-helper --normal-tumor-bam-pair \"$line\" --roi-file \"$roi_file\" ".
    "--reference-sequence \"$ref_seq\" --output-file \"$roi_covg_dir/$sample.covg\"$optional_params\'";

    # If user only wants the calcRoiCovg commands, write them to file and skip running calcRoiCovg
    if( defined $cmd_list_file )
    {
      $calcRoiCovg_cmd = $cmd_prefix . " $calcRoiCovg_cmd" if( defined $cmd_prefix );
      $cmdFh->print( "$calcRoiCovg_cmd\n" );
      next;
    }

    # If the calcRoiCovg output was already generated, then don't rerun it
    if( -s "$roi_covg_dir/$sample.covg" )
    {
      print "$sample.covg found in $roi_covg_dir. Skipping re-calculation.\n";
    }
    # Run the calcRoiCovg command on this tumor-normal pair. This could take a while
    else {
      my %params = (
        normal_tumor_bam_pair => $line,
        roi_file => $roi_file,
        reference_sequence => $ref_seq, 
        output_file => $roi_covg_dir."/".$sample.".covg",
      );
      if ($normal_min_depth) {
        $params{"normal_min_depth"} = $normal_min_depth;
      }
      if ($tumor_min_depth) {
        $params{"tumor_min_depth"} = $tumor_min_depth;
      }
      if ($min_mapq) {
        $params{"min_mapq"} = $min_mapq;
      }
      my $cmd = Genome::Model::Tools::Music::Bmr::CalcCovgHelper->create(%params);
      my $rv = $cmd->execute;
      if(!$rv)
      {
        print STDERR "Failed to execute: $calcRoiCovg_cmd\n";
        next;
      }
      else
      {
        print "$sample.covg generated and stored to $roi_covg_dir.\n";
      }
    }

    # Read the calcRoiCovg output and count covered bases per gene
    my %geneCovg = ();
    my ( $tot_covd, $tot_at_covd, $tot_cg_covg, $tot_cpg_covd );
    my $roiCovgFh = IO::File->new( "$roi_covg_dir/$sample.covg" );
    while( my $line = $roiCovgFh->getline )
    {
      chomp( $line );
      if( $line =~ m/^#NonOverlappingTotals/ )
      {
        ( undef, undef, undef, $tot_covd, $tot_at_covd, $tot_cg_covg, $tot_cpg_covd ) = split( /\t/, $line );
      }
      elsif( $line !~ m/^#/ )
      {
        my ( $gene, undef, $length, $covd, $at_covd, $cg_covd, $cpg_covd ) = split( /\t/, $line );
        $geneCovg{$gene}{len} += $length;
        $geneCovg{$gene}{covd_len} += $covd;
        $geneCovg{$gene}{at} += $at_covd;
        $geneCovg{$gene}{cg} += $cg_covd;
        $geneCovg{$gene}{cpg} += $cpg_covd;
      }
    }
    $roiCovgFh->close;

    # Write the per-gene coverages to a file named after this sample_name
    my $geneCovgFh = IO::File->new( "$gene_covg_dir/$sample.covg", ">" );
    $geneCovgFh->print( "#Gene\tLength\tCovered\tAT_covd\tCG_covd\tCpG_covd\n" );
    foreach my $gene ( sort keys %geneCovg )
    {
      $geneCovgFh->print( join( "\t", $gene, $geneCovg{$gene}{len}, $geneCovg{$gene}{covd_len},
         $geneCovg{$gene}{at}, $geneCovg{$gene}{cg}, $geneCovg{$gene}{cpg} ), "\n" );
    }
    $geneCovgFh->close;

    # Write total coverages for this sample to a file
    $totCovgFh->print( "$sample\t$tot_covd\t$tot_at_covd\t$tot_cg_covg\t$tot_cpg_covd\n" );
  }
  $bamFh->close;
  $cmdFh->close if( defined $cmd_list_file );
  $totCovgFh->close unless( defined $cmd_list_file );

  return 1;
}

1;
