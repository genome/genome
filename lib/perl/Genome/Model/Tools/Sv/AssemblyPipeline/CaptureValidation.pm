package Genome::Model::Tools::Sv::AssemblyPipeline::CaptureValidation;

use strict;
use warnings;
use Carp;
use Genome;
use Statistics::Descriptive;

class Genome::Model::Tools::Sv::AssemblyPipeline::CaptureValidation {
  is => 'Command',
  has => [
    assembly_file => { is => 'Text', doc => "Assembled FASTA sequences of the SVs to validate" },
    sv_file => { is => 'Text', doc => "Describes an SV per line. Use BreakDancer output format" },
    normal_bam => { is => 'Text', doc => "Path to the Normal BAM used when calling SVs" },
    tumor_bams => { is => 'Text', doc => "Comma-separated tissue names and paths to one or more tumor BAMs. E.g. tumor:path,mets:path,xeno:path" },
    output_file => { is => 'Text', doc => "Path to the resulting annotated SV output" },
    dump_read_alignments => { is => 'Boolean', doc => "Enable to see the read alignments for SV-specific reads", is_optional => 1, default => 0 },
    no_duplicate_read_sequence => { is => 'Boolean', doc => "Enable to skip printing duplicate sequences", is_optional => 1, default => 1 },
    extend => { is => 'Integer', doc => "Number of bases beyond breakpoint that the read has to span", is_optional => 1, default => 10 },
    max_unaligned_end_bps => { is => 'Integer', doc => "Max allowed unaligned bases at end of a read", is_optional => 1, default => 1 },
    buffer => { is => 'Integer', doc => "Flank (bps) around breakpoints to fetch reads from; 2*buffer fetched for refseq", is_optional => 1, default => 500 },
    min_fraction_diff => { is => 'Number', doc => "Min fraction of difference between alignments to SV contig vs refseq", is_optional => 1, default => 0.34 },
    max_percent_subs => { is => 'Number', doc => "", is_optional => 1, default => 1 },
    max_percent_indels => { is => 'Number', doc => "", is_optional => 1, default => 1 },
    min_score => { is => 'Integer', doc => "", is_optional => 1, default => 50 },
    percent_contamination => { is => 'Number', doc => "Currenly unused", is_optional => 1, default => 0 },
  ],
};

sub help_brief {
  "Given BreakDancer output, produces files for capture validation";
}

# Example:
# bsub -R 'select[mem>8000] rusage[mem=8000]' gmt sv assembly-pipeline capture-validation --sv-file /gscmnt/sata149/info/medseq/AML/Ley-M1-M3-Capture-Validation/AML10/AML10.allchr.Q40.somatic.assembled.HQfiltered.csv --assembly-file /gscmnt/sata149/info/medseq/AML/Ley-M1-M3-Capture-Validation/AML10/AML10.allchr.Q40.somatic.assembled.HQfiltered.fasta --normal-bam /gscmnt/sata877/info/model_data/2859174445/build103345592/alignments/103345592_merged_rmdup.bam --tumor-bams tumor:/gscmnt/sata878/info/model_data/2859174519/build103345660/alignments/103345660_merged_rmdup.bam --output-file /gscmnt/sata149/info/medseq/AML/Ley-M1-M3-Capture-Validation/AML10/allchr.110404.csv

sub help_detail {
  return <<HELP;
Requires BreakDancer file, assembly file, and bam files as input. Resulting output file contains
BreakDancer lines along with the following additional columns:
Fasta header of SV contig, and for each tissue:
- total number of reads mapping to +/- 500 bp of each breakpoint on reference
- total number of SV-supporting reads

Execution might require up to 8GB of memory. Submit to LSF as follows:
bsub -R 'select[mem>8000] rusage[mem=8000]' gmt sv assembly-pipeline capture-validation ...
HELP
}

sub execute {
  my $self = shift;

  # Grab arguments
  my $assemblyFastaFile = $self->assembly_file;
  my $svFile = $self->sv_file;
  my $normalBam = $self->normal_bam;
  my $tumorBams = $self->tumor_bams;
  my $outputFile = $self->output_file;
  my $dumpReadAlignments = $self->dump_read_alignments;
  my $noDuplicateReadSequence = $self->no_duplicate_read_sequence;
  my $extend = $self->extend;
  my $maxUnalignedEndBps = $self->max_unaligned_end_bps;
  my $maxPercentSubs = $self->max_percent_subs;
  my $maxPercentIndels = $self->max_percent_indels;
  my $minScore = $self->min_score;
  my $percentContamination = $self->percent_contamination;
  my $buffer = $self->buffer;
  my $minFractionDiff = $self->min_fraction_diff;

  # Use -masklevel 101 to get all hits, but ReadRemap::createHitObjects only returns one hit per query
  # Ken's parameters do not give a complete alignment. Seems to penalize too much for mismatches at ends of reads
  # Ken's parameters: -bandwidth 20 -minmatch 20 -minscore 25 -penalty -10 -discrep_lists -tags -gap_init -10 -gap_ext -1
  # my $crossMatchParameters = "-bandwidth 20 -minmatch 20 -minscore 25 -penalty -10 -tags -gap_init -10 -gap_ext -1  -masklevel 101";
  my $crossMatchParameters = "-discrep_lists -minmatch 10 -maxmatch 10 -minscore 15";

  # Parse out the paths to the BAMs and store them in a hash using tissue names as keys
  my %tissueToBamFile = map { split( /\:/ ) } split( /,/, $tumorBams );
  $tissueToBamFile{normal} = $normalBam;

  # Check parameters.
  my $error = "";
  $error .= "\n  Assembly fasta file was not found" unless( -s $assemblyFastaFile );
  $error .= "\n  Sv file was not found" unless( -s $svFile );
  foreach( keys %tissueToBamFile ) {
    $error .= "\n  BAM file for $_ tissue not found" unless( -s $tissueToBamFile{$_} );
  }
  ( $error eq "" ) or die "Halted execution due to following:$error\n";

  # Parse SV file to get hash of regions and hash of fasta header IDs
  my ( @entireSvFile, $line, $chrA, $chrB, $bpA, $bpB, $id, %regions, %ids );
  open( IN, "<$svFile" ) or die "Could not open $svFile: $!";
  @entireSvFile = <IN>;
  chomp( @entireSvFile );
  close( IN );
  foreach $line( @entireSvFile ) {
    next if( $line =~ /^#/ );

    #ID     CHR1    OUTER_START     INNER_START     CHR2    INNER_END       OUTER_END       TYPE    ORIENTATION     MINSIZE MAXSIZE SOURCE  SCORES  Copy_Number
    #9.1     9       21911178        21911178        9       21971899        21971899        DEL     +-      60720   60720   TEST    402     NA      NA      NA

    #parse merged bd file line
    my ($id,$chrA,$bpA,undef,$chrB,$bpB,undef,$type,$orientation,$minsize,$maxsize,$source,$score) = split /\t/,$line;

#    my $bdRef = Genome::Model::Tools::Sv::AssemblyPipeline::BreakDancerLine::newline( $self,$line );
#    ( $chrA, $bpA, $chrB, $bpB ) = $bdRef->chromosomesAndBreakpoints();
    ( defined $bpA && $bpA =~ /^\d+$/ && defined $bpB && $bpB =~ /^\d+$/ ) or
      confess "Did not get chr and breakpoints from '$line'";
#    $id = $bdRef->Id();
    $regions{"$chrA.$bpA.$chrB.$bpB"} = 1;
    $ids{$id} = 0;
  }
  my ( $regionsRef, $idRef ) = ( \%regions, \%ids );

  # Temp file names to keep SV contig and reference sequences
  my $contigSequenceFile = "/tmp/contigSequenceFile".rand();
  my $refSequenceFile = "/tmp/refSequenceFile".rand();

  # Get SV contig sequences and put into *fasta file
  # This just goes through and pulls out the contigs that are listed in SV file.  It should be a
  # subset of the total contigs found in $assemblyFastaFile.  The point is to only compare reads
  # to the assembly contigs of interest. At the same time, the fasta header to ID look-up hash is made
  $idRef = Genome::Model::Tools::Sv::AssemblyPipeline::ReadRemap::getAssemblySequences( $idRef, $assemblyFastaFile, $contigSequenceFile );

  # Get regions surrounding each SV breakpoint from Build36 reference and put into a fasta file
  Genome::Model::Tools::Sv::AssemblyPipeline::ReadRemap::getBuild36ReferenceSequences( $regionsRef, $refSequenceFile, 2*$buffer );

  # Make sure the sequences files exist with non-zero size
  ( -s $contigSequenceFile && -s $refSequenceFile ) or die "Did not get contig sequence and/or reference sequence file";

  # Go through each entry of SV file and find reads that cross breakpoint in SV contig
  open( OUT, ">$outputFile" ) or die "Could not open '$outputFile': $!";
  foreach my $line ( @entireSvFile ) {
    print OUT "$line";
    if( $line =~ /^#/ ) {
      print OUT "\n";
      next;
    }

    my ($id,$chrA,$bpA,undef,$chrB,$bpB,undef,$type,$orientation,$minsize,$maxsize,$source,$score) = split /\t/,$line;
    #my $bdRef = Genome::Model::Tools::Sv::AssemblyPipeline::BreakDancerLine::newline($self,$line );
    #my $id = $bdRef->Id();
    my $fastaHeader = $$idRef{$id};
    $fastaHeader =~ s/\>//;
    print OUT "\t$fastaHeader";
    if( $fastaHeader !~ /Var\:(\w+)\.(\d+)\.(\w+)\.(\d+)\.(\w+)\./ ||
        $fastaHeader !~ /Ins\:(\d+)\-(\-|\d+)/ ) {
      print OUT "no fasta sequence\n";
      next;
    }

    foreach my $tissue ( sort keys %tissueToBamFile ) {
      my $bamFile = $tissueToBamFile{$tissue};

      # All reads aligned to assembly contigs and reference sequence
      my ( $readCount, $uniqueReadCount, $hitsToAssemblyRef, $hitsToNormalRef ) =
        Genome::Model::Tools::Sv::AssemblyPipeline::ReadRemap::remapByCrossMatch( $line, $buffer, $fastaHeader, $bamFile, $contigSequenceFile,
                                      $refSequenceFile, $crossMatchParameters, $noDuplicateReadSequence );

      # If there are no capture reads from the region, can't do any more
      if ( $readCount == 0 ) { print OUT "\t$tissue.totalReads:0\t$tissue.SvReadCount:0"; next; }

      # Get SV-specific read alignments based on cross_match realignment
      my ( undef, $crossesSvBreakpointRef ) =
        Genome::Model::Tools::Sv::AssemblyPipeline::ReadRemap::svSpecificHits( $fastaHeader, $hitsToAssemblyRef, $hitsToNormalRef, $extend,
                                   $maxUnalignedEndBps, $maxPercentSubs, $maxPercentIndels, $minScore );
      my $filteredSvReadRef =
        Genome::Model::Tools::Sv::AssemblyPipeline::ReadRemap::removeMarginalSvHits( $crossesSvBreakpointRef, $hitsToNormalRef, $minFractionDiff );
      print OUT "\t$tissue.totalReads:$readCount\t$tissue.SvReadCount:", scalar( keys %{$filteredSvReadRef} );

      # For read alignment dump
      if( $dumpReadAlignments ) {
        print OUT "\n$tissue\n";
        Genome::Model::Tools::Sv::AssemblyPipeline::ReadRemap::displayAlignments( $filteredSvReadRef, $hitsToNormalRef );
      }

      # When the SAM reads are converted to *.fasta, they are not complemented based on flag
      # so this is bogus.....
      my $fractionComplemented = Genome::Model::Tools::Sv::AssemblyPipeline::ReadRemap::fractionComplemented( $filteredSvReadRef );
      print OUT "\t$tissue.Complemented.$fractionComplemented";
    }
    print OUT "\n";
  }

  close( OUT );
  unlink $contigSequenceFile;
  unlink $refSequenceFile;

  return 1;
}

1;
