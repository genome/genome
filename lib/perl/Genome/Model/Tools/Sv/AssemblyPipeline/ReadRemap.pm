package Genome::Model::Tools::Sv::AssemblyPipeline::ReadRemap;

class Genome::Model::Tools::Sv::AssemblyPipeline::ReadRemap {
        is => 'Genome::Model::Tools::Sv::AssemblyPipeline'
};

#
#   Run on 64-bit machine
#
#
#
#
#

use strict;
use warnings;

use Carp;
use FindBin qw($Bin);
use lib "$FindBin::Bin";
use PostData;
#use lib "/gscuser/jwallis/svn/perl_modules/test_project/jwallis";
#use Genome::Model::Tools::Sv::AssemblyPipeline::Hits;

# This is used to parse the cross_match hits
my $number = "\\d+\\.?\\d*";
my $deleted = "\\(\\s*\\d+\\s*\\)";
my $AlignmentLine = "($number)\\s+($number)\\s+($number)\\s+($number)\\s+(\\S+)\\s+(\\d+)\\s+(\\d+)\\s+($deleted)\\s+(C\\s+)?(\\S+)\\s+($deleted)?\\s*(\\d+)\\s+(\\d+)\\s*($deleted)?";



sub getReads {
    # Input: '$chr' (fasta header used to create *.bam file)
    #         $start, $stop (reads cross ANY bases between $start, $stop inclusive)
    #         $bamFile
    #         $noDuplicates -- set to 1 to get de-duplicated reads; 
    # Return: ref to array of reads that map to anywhere defined by $chr, $start, $stop OR mates of reads that do
    #
    # No minimum quality value is used so mates of mapped reads can be captured.
    # The following gets mapped and unmapped reads that are not marked as duplicates:
    #       samtools view -F 0x0400 <bam> <chr:start-stop>
    #       (This only gets reads that are mapped: samtools view -F 0x0404 <bam> <chr:start-stop> -- don't want this for SVs)

    my ($chr, $start, $stop, $bamFile, $noDuplicates) = @_;
    (-e $bamFile) || confess "Could not find '$bamFile': $!";

    # See if duplicate reads are allowed
    if ( !defined $noDuplicates ) { $noDuplicates = 0; }
    my $duplicateFlag = "";
    if ( $noDuplicates ) { $duplicateFlag = "-F 0x0400"; }

    open(SAM, "samtools view $duplicateFlag $bamFile $chr:$start-$stop |") || confess "Could not open pipe to samtools: $!";
    my @breakPointReads = <SAM>;
    close SAM;
    return \@breakPointReads;
}

sub getAssemblySequences {
    # Input: ref to hash of unique ID(s), fasta sequence file, output file
    # Write given sequences to output file
    # Update the values for unique ID hash to get the fasta header

    my ($idRef, $fastaFile, $outputFile ) = @_;
    my ( $reading, @entireFile, $line, $sequence, $sequencesFound, );
    $reading = 0;
    open(IN, "< $fastaFile") || confess "Could not open '$fastaFile': $!";
    @entireFile = <IN>;
    open(OUT, "> $outputFile") || confess "Could not open '$outputFile' for output: $!";
    foreach $line (@entireFile) {
        chomp $line;
        if ( $line =~ />(\w+\.\d+)\,/ ) {
            if ( defined $$idRef{$1} ) { $reading = 1; $$idRef{$1} = $line; } else { $reading = 0; }
        }
        if ( $reading ) { print OUT "$line\n"; }
    }
    return $idRef;
}

sub getBuild36ReferenceSequences {
    # Input: ref to hash with key = "$chrA.$bpA.$chrB.$bpB", output file, $buffer
    # Use 'expiece' to write region +/- $buffer of each position
    # Chromosome-specific reference sequences are at:
    # /gscmnt/sata180/info/medseq/biodb/shared/Hs_build36/Homo_sapiens.NCBI36.45.dna.chromosome.$chr.fa

    my ($positionRef, $outputFile, $buffer) = @_;
    my ( $position, $remove, $chrA, $bpA, $chrB, $bpB, $start, $stop, @output, $seq, );
    if (!defined $buffer) { $buffer = 500; }

    open(OUT, "> $outputFile") || confess "Could not open '$outputFile' for output: $!";

    # This is what the fasta header looks like when using 'expiece'.
    # Want to remove the full path to the chromosome
    $remove = "\/gscmnt\/sata180\/info\/medseq\/biodb\/shared\/Hs_build36\/Homo_sapiens.NCBI36.45.dna.chromosome.";

    foreach $position ( sort keys %{$positionRef} ) {
        if ( $position =~ /(\w+)\.(\d+)\.(\w+)\.(\d+)/ ) {
            $chrA = $1; $bpA = $2; $chrB = $3; $bpB = $4;
        } else {
            confess "Unexpected format for position: '$position'";
        }

        # Do first coordinate
        $start = $bpA - $buffer;
        $stop = $bpA + $buffer;
        my $chrFile = "/gscmnt/sata180/info/medseq/biodb/shared/Hs_build36/Homo_sapiens.NCBI36.45.dna.chromosome.$chrA.fa";
        open(EXP, "expiece $start $stop $chrFile |") || confess "Could not open pipe for 'expiece $start $stop $chrFile'";
        @output = <EXP>;
        close EXP;
        foreach my $seq (@output) {
            chomp $seq;
            if ( $seq =~ />/ ) { 
                $seq =~ s/$remove//;
                $seq =~ s/fa from $start to $stop/$start.$stop/;
            }
            print OUT "$seq\n";
        }

        # Do second coordinate
        $start = $bpB - $buffer;
        $stop = $bpB + $buffer;
        $chrFile = "/gscmnt/sata180/info/medseq/biodb/shared/Hs_build36/Homo_sapiens.NCBI36.45.dna.chromosome.$chrB.fa";
        open(EXP, "expiece $start $stop $chrFile |") || confess "Could not open pipe for 'expiece $start $stop $chrFile'";
        @output = <EXP>;
        close EXP;
        foreach my $seq (@output) {
            chomp $seq;
            if ( $seq =~ />/ ) { 
                $seq =~ s/$remove//;
                $seq =~ s/fa from $start to $stop/$start.$stop/;
            }
            print OUT "$seq\n";
        }


    } # matches 'foreach $position'
    close OUT;
}
sub getBuild37ReferenceSequences {
    # Input: ref to hash with key = "$chrA.$bpA.$chrB.$bpB", output file, $buffer
    # Use 'expiece' to write region +/- $buffer of each position
    # Chromosome-specific reference sequences are at:
    # /gscmnt/ams1102/info/model_data/2869585698/build106942997/$chr.fa

    my ($positionRef, $outputFile, $buffer) = @_;
    my ( $position, $remove, $chrA, $bpA, $chrB, $bpB, $start, $stop, @output, $seq, );
    if (!defined $buffer) { $buffer = 500; }

    open(OUT, "> $outputFile") || confess "Could not open '$outputFile' for output: $!";

    # This is what the fasta header looks like when using 'expiece'.
    # Want to remove the full path to the chromosome
    $remove = "\/gscmnt\/ams1102\/info\/model_data\/2869585698\/build106942997\/";
    #$remove = "\/gscmnt\/sata180\/info\/medseq\/biodb\/shared\/Hs_build36\/Homo_sapiens.NCBI36.45.dna.chromosome.";

    foreach $position ( sort keys %{$positionRef} ) {
        if ( $position =~ /(\w+)\.(\d+)\.(\w+)\.(\d+)/ ) {
            $chrA = $1; $bpA = $2; $chrB = $3; $bpB = $4;
        } else {
            confess "Unexpected format for position: '$position'";
        }

        # Do first coordinate
        $start = $bpA - $buffer;
        $stop = $bpA + $buffer;
        my $chrFile = "/gscmnt/ams1102/info/model_data/2869585698/build106942997/$chrA.fa";
        open(EXP, "expiece $start $stop $chrFile |") || confess "Could not open pipe for 'expiece $start $stop $chrFile'";
        @output = <EXP>;
        close EXP;
        foreach my $seq (@output) {
            chomp $seq;
            if ( $seq =~ />/ ) { 
                $seq =~ s/$remove//;
                $seq =~ s/fa from $start to $stop/$start.$stop/;
            }
            print OUT "$seq\n";
        }

        # Do second coordinate
        $start = $bpB - $buffer;
        $stop = $bpB + $buffer;
        $chrFile = "/gscmnt/ams1102/info/model_data/2869585698/build106942997/$chrB.fa";
        open(EXP, "expiece $start $stop $chrFile |") || confess "Could not open pipe for 'expiece $start $stop $chrFile'";
        @output = <EXP>;
        close EXP;
        foreach my $seq (@output) {
            chomp $seq;
            if ( $seq =~ />/ ) { 
                $seq =~ s/$remove//;
                $seq =~ s/fa from $start to $stop/$start.$stop/;
            }
            print OUT "$seq\n";
        }


    } # matches 'foreach $position'
    close OUT;
}


sub samReadPassesFilter {
    # Input: $read in SAM format
    #        $maxMismatch
    #        $maxSoftClip 
    # Return: 1 if read has <= $maxMismatch and <= $maxSoftClip bases

    my ($read, $maxMismatch, $maxSoftClip) = @_;
    (defined $read && defined $maxMismatch && defined $maxSoftClip) || carp
        "Input parameter(s) not defined: '$read',  '$maxMismatch', '$maxSoftClip'";

    my $passesFilter = 1;

    # Mismatches
    if ( $read =~ /NM:i:(\d+)/ && $1 > $maxMismatch ) { $passesFilter = 0; }   
    # Softclip
    my (undef, undef, undef, undef, undef, $cigar) = split /\s+/, $read;
    if ( $cigar =~ /(\d+)S/ && $1 > $maxSoftClip ) { $passesFilter = 0; }
    if ( $cigar =~ /(\d+)S\w+\D(\d+)S/ && $1+$2 > $maxSoftClip ) { $passesFilter = 0; }

    return $passesFilter;
}

sub samReadCrossesBreakpoints {
    # Input: aligned $read in SAM format
    #        $start, $stop  -- region the read should cover (expressed in coordinates of $chr to which the read was aligned)
    # Return: 1 if read covers both coordinates.  It can't stop at a coordinate; it must extend by at least one base
    # The read must be mapped, not just it's mate

    my ($read, $start, $stop) = @_;
    (defined $read && defined $start && defined $stop) || carp
        "Input parameter(s) not defined: '$read', '$start', '$stop'";

    my ( $flag, $readStart, $readStop, $cigar, $seq );
    (undef, $flag, undef, $readStart, undef, $cigar, undef, undef, undef, $seq) = split /\s+/, $read;

    # Read should be mapped  0x0004 => the query sequence itself is unmapped
    if ( $flag & 0x0004) { return 0; }

    $readStop = $readStart + length($seq) - 1;

    return ( $readStart < $start && $readStop > $stop )
}

sub convertSamToFasta {
    # Input: ref to array of reads in SAM format
    #        File name for output
    #        Flag set to 1 if a quality file 'outFile.qual' should be written
    # Output: Fasta sequence where read names are "$name.$cigar.$flag"
    #         Quality file if flag set to 1
    # Return: number of sequences written


    my ( $readRef, $outFile, $writeQualityFile ) = @_;
    if ( !defined $writeQualityFile ) { $writeQualityFile = 0; }

    my ( $line, $name, $flag, $cigar, $seq, %allReadNames, $qualityString,  );
    %allReadNames = ();
    open(FASTA, "> $outFile") || confess "Could not open '$outFile': $!";
    if ( $writeQualityFile ) { open(QUAL, "> $outFile.qual") || confess "Could not open '$outFile.qual': $!"; }
    foreach $line ( @{$readRef} ) {
        chomp $line;
        ($name, $flag, undef, undef, undef, $cigar, undef, undef, undef, $seq, $qualityString) = split /\s+/, $line;

        # Make $name unique 
        $name = "$name.$cigar.$flag";
        ( !defined $allReadNames{$name} ) ||
            confess "Duplicate read name '$name'";
        $allReadNames{$name} = 1;
        print FASTA ">$name\n$seq\n";
        if ( $writeQualityFile ) {
            $qualityString = fastqToPhred($qualityString);
            print QUAL ">$name\n$qualityString\n";
        }
    }

    return ( scalar(keys %allReadNames) );
}

sub fastqToPhred {
    # Convert fastq quality values to phred quality values
    my $qualityLine = $_[0];
    my $phredLine = "";
    foreach my $q (split //, $qualityLine) {
        my $quality = ord($q) - 33;
        $phredLine .= "$quality ";
    }
    return $phredLine;
}

sub runCrossMatch {
    # Input: query file, subject file, parameters (as string), optional output file
    # Return: ref to array of cross_match results.  Will also write to file if file is given

    my ( $query, $subject, $parameters, $outFile ) = @_;
    (-e $query) || confess "'$query' does not exist"; 
    (-e $subject) || confess "'$subject' does not exist";
    (defined $parameters) || confess "cross_match parameters are not defined";
    if ( !defined $outFile ) { $outFile = "/dev/null"; }


    open(CROSS, "cross_match $parameters $query $subject |") || 
        confess "Not able to start 'cross_match $parameters $query $subject'";
    my @output = <CROSS>;
    close CROSS;
    open(OUT, "> $outFile") || confess "Cold not open '$outFile': $!";
    print OUT "@output";
    close OUT;

    return \@output;
}


sub createHitObjects {
    # Input: either a file with cross_match results or ref to array of cross_match results
    # Return: ref to hash with key = query name; value = Hits object
    #   NOTE: only one hit is returned for each query.  It is the hit with the best score.

    my $input = $_[0];
    my ( @allCrossMatch, $line, $hit, $query, %hitList, );

    # Get array with all cross_match output
    if ( -e $input ) {
        open(IN, "< $input") || confess "Could not open '$input': $!";
        @allCrossMatch = <IN>;
        close IN;
    } else {
        @allCrossMatch = @{$input};
    }

    foreach $line ( @allCrossMatch ) {
        chomp $line;
        if ( $line =~ /$AlignmentLine/ ) {
            $hit = Genome::Model::Tools::Sv::AssemblyPipeline::Hits->new;
            $hit->addCrossMatchLine($line);
            $query = $hit->queryName();
            # It this query already has a hit, choose the one with the highest score
            if ( !defined $hitList{$query} || $hit->score() > $hitList{$query}->score() ) {
                $hitList{$query} = $hit;
            }
        }

        # Output from -discrep_lists follows the alignment line for a given hit
        # I     847  T(40)   1260  cctcagTcctggg
        # I-5   913  CATCC(40)   1290  ccatccCATCCacccac
        # D-24   586  T(40)   1040  tgatttTtatata  
        # D-60   651  T(40)   1107  gatggaTagatgg
        # D-$size  $queryStart T(40)  $subjectStart
        my ( $indel, $length, $queryStart, $targetStart );
        if ( $line =~ /\s*(I|D)\-?(\d*)\s+(\d+)\s+\S+\s+(\d+)/ ) {
            $indel = $1; $length = $2; $queryStart = $3; $targetStart = $4;
            if ( !defined $length || $length eq "" ) { $length = 1; } 
            if ( $indel eq "I" ) {
                $hitList{$query}->addInsertion($queryStart, $targetStart, $length);
            } elsif ( $indel eq "D" ) {
                $hitList{$query} ->addDeletion($queryStart, $targetStart, $length);
            }
        }
    }

    return \%hitList;
}
 
sub crossMatchHitPassesFilter {
    # Input: Hits.pm object
    # Return: 1 if passes criteria

    my ($hitObj, $maxUnalignedEndBases, $maxPercentSubs, $maxPercentIndels, $minScore) = @_;
    (defined  $maxUnalignedEndBases && defined $maxPercentSubs && defined $maxPercentIndels && defined $minScore) ||
        confess "input parameters not defined";
    if ( $hitObj->pastQueryEnd() > $maxUnalignedEndBases || 
         $hitObj->queryStart() > $maxUnalignedEndBases ||
         $hitObj->percentSubs() > $maxPercentSubs || 
         $hitObj->score() < $minScore || 
         $hitObj->percentDels() > $maxPercentIndels ||
         $hitObj->percentInserts() > $maxPercentIndels 
        ) { return 0; }

    return 1;
}

sub crossMatchHitCrossesBreakpoints {
    # Input: Hits.pm object, $start, $stop
    # Return: 1 if read covers both coordinates.  It can't stop at a coordinate; it must extend by at least one base

    ####  
    #  For the case where the read stops a few bases past a coordinate, the unaligned portion is relevant (pastQueryEnd and 
    #  queryStart)


    my ($hitObj, $start, $stop) = @_;
    # Make sure $start < $stop
    if ( $start > $stop ) { ($start, $stop) = ($stop, $start); }
    (defined $start && defined $stop && $start =~ /^\-?\d+$/ && $stop =~ /^\d+$/) ||
        confess "Expected integer start and stop";

    my $subjectStart = $hitObj->subjectStart();
    my $subjectEnd = $hitObj->subjectEnd();
    if ( $subjectStart > $subjectEnd ) { ($subjectStart, $subjectEnd) = ($subjectEnd, $subjectStart); }

    return ( $subjectStart < $start && $subjectEnd > $stop );

}

sub removeMarginalSvHits {
    # Input: Ref to hash of Hit objects of reads to SV contig
    #        ref to hash of Hit objects of reads to reference sequence
    #        maximum fraction difference
    #    ** The keys for the two hashes must be same so we are comparing alignments of a single read to SV contig and reference **
    # Return: Ref to hash of Hit objects with SV-supporting hits that are not close to 
    #         their alignments to reference.
    # The purpose of this is to remove SV-supporting reads that just barely pass the filters
    # e.g. number of mismatches for read to reference could be 5 and number of mismatches of read
    # to SV contig could be 4.

    my ( $hitsToSvRef, $hitsToReferenceRef, $minFractionDifference ) = @_;
    my ( %filteredHits, $readName, );
    foreach $readName (keys %{$hitsToSvRef}) {
        if ( defined $$hitsToReferenceRef{$readName} && 
             &crossMatchHitsSimilar( $$hitsToSvRef{$readName},$$hitsToReferenceRef{$readName},$minFractionDifference)
            ) {
            # Don't use this read, it's alignment to the SV contig is too close to alignment to reference
            next;
        }

        # This read is OK.
        $filteredHits{$readName} = $$hitsToSvRef{$readName};
    }

    return \%filteredHits;
}



sub crossMatchHitsSimilar {
    # Input: Two hit objects and min fraction difference
    # Return: 0 or 1 depending if %subs, %dels, %inserts or unaligned bases at ends
    #         differ by more than max fraction difference
    # Used to compare alignment of read to reference and assembly contig
    #    
    my ( $hit_1, $hit_2, $fractionDiff ) = @_;

    return ( fractionDifference($hit_1->pastQueryEnd(), $hit_2->pastQueryEnd()) < $fractionDiff  &&
             fractionDifference($hit_1->queryStart(), $hit_2->queryStart()) < $fractionDiff &&
             fractionDifference($hit_1->percentSubs(),  $hit_2->percentSubs()) < $fractionDiff &&
             fractionDifference($hit_1->percentDels(),  $hit_2->percentDels()) < $fractionDiff &&
             fractionDifference($hit_1->percentInserts(),  $hit_2->percentInserts()) < $fractionDiff
        );   
}

sub fractionDifference {
    # Return fraction difference between two numbers
    # If one of the numbers is 0, add 0.1 so fraction difference isn't huge

    my ( $num1, $num2 ) = @_;
    if ($num1 == $num2) { return 0; }

    if ($num1 > $num2) { ($num1, $num2) = ($num2, $num1); }

    if ( $num1 == 0 ) { $num1 = 0.1; }
    return ( ($num2 - $num1)/$num2 );
}


sub uniqueCrossMatchAlignments {
    # Input: Ref to hash of Hits objects
    # Return: Ref to hash of Hits objects that have unique alignments.  If alignments are the
    #         same, one is arbitrarily chosen to represent all
    #         Key = query name, value = ref to Hits object
    # A hash is make where key is based on alignment line without the query name, value = alignment line
    # Then Hits objects are made from the unique alignments

    my $hitObjRef = $_[0];
    my (%newHitObjects, $line, @cols, $hashKey, %uniqueAlignments, $hit, );

    # This makes hash with key = alignment line without read name
    foreach my $id (keys %{$hitObjRef}) {
        $hashKey = "";
        $line = $$hitObjRef{$id}->alignmentLine();
        $line =~ s/^\s*//; $line =~ s/\s*$//;
        @cols = split /\s+/, $line;

        # The query name is in the 5th column ($i = 4)
        for (my $i = 0; $i <= $#cols; $i++ ) {
            if ( $i != 4 ) { $hashKey .= "$cols[$i]"; }
        }
        $uniqueAlignments{$hashKey} = $line;
    }

    # Now make a new hash of Hits objects
    foreach my $id ( keys %uniqueAlignments ) {
        $line = $uniqueAlignments{$id};
        ($line =~ /$AlignmentLine/ ) || confess "'$line' is not an alignment line";
        $hit = Genome::Model::Tools::Sv::AssemblyPipeline::Hits->new;
        $hit->addCrossMatchLine($line);
        my $query = $hit->queryName();
        $newHitObjects{$query} = $hit;
    }

    return \%newHitObjects;
}


sub uniqueSamSequenceReads {
    # Input: ref to array of reads in SAM format
    # Return: ref to array of reads that have unique sequence
    # 
    #    NOTE: this ignores pairs so is not a good way to de-duplicate
    #    if aligner takes into account paired reads
    #    Essentially turns reads into fragment reads.  Just uses hash to 
    #    choose the last one of the unique sequences

    my $readRef = $_[0];
    my ( $line, %uniqueReads, @unique, $seq );
    foreach $line ( @{$readRef} ) {
        chomp $line;
        (undef, undef, undef, undef, undef, undef, undef, undef, undef, $seq) = split /\s+/, $line;
        $uniqueReads{$seq} = $line;
    }
    @unique = keys %uniqueReads;
    return \@unique;
}

sub isCrossMatchAlignmentLine {
    # Returns 0 or 1 depending on whether line matches expected cross_match ouput line
    my $line = $_[0];
    if ( !defined $line || $line !~ /\w+/ ) { return 0; }
    return ( $line =~ /$AlignmentLine/ );
}


sub  remapByCrossMatch {
    # Input: ( $line, $buffer, $fastaHeader, $bamFile, $contigSequenceFile, $referenceSequenceFile, $crossMatchParameters, $noDuplicates )
    #       BreakDancer line, +/-$buffer around breakpoints from which to get reads
    #       $noDuplicates: set to 1 if reads should be de-duplicated; default is 0
    #
    # Align all reads to SV contig and reference.  
    #
    # Return: ($readCount, $uniqueReadCount, $hitsToAssemblyRef, $hitsToNormalRef)
    #     number of reads and mates that map +/- buffer to SV breakpoints defined in BreakDancer line
    #     number of above reads that are unique
    #     All hits that align to SV assembly contig(s): ref to hash with key = read name; value = Hits object 
    #     All hits that align to reference sequence: ref to hash with key = read name; value = Hits object 
    #

    my ( $line, $buffer,  $fastaHeader, $bamFile, $contigSequenceFile, $referenceSequenceFile, $crossMatchParameters, $noDuplicates ) = @_;

    my ($id,$chrA,$bpA,undef,$chrB,$bpB,undef,$type,$orientation,$minsize,$maxsize,$source,$score) = split /\t/,$line;
    #my $bdRef = BreakDancerLine->new($line);
    #my ($chrA, $bpA, $chrB, $bpB) = $bdRef->chromosomesAndBreakpoints();
    defined ( $bpA && $bpA =~ /^\d+$/ && defined $bpB && $bpB =~ /^\d+$/ ) ||
        confess "Did not get chr and breakpoints from '$line'";

    my ( $readCount, $uniqueReadCount, $crossMatchResults, $hitsToAssemblyRef, $hitsToNormalRef );

    # Get all reads surrounding each breakpoint. The regions may overlap with buffer so make sure reads are unique
    my $readRef = &getReads($chrA, $bpA-$buffer, $bpA+$buffer, $bamFile, $noDuplicates);
    my %uniqueEntries = ();
    foreach ( @{$readRef} ) { $uniqueEntries{$_} = 1; }
    $readRef = &getReads($chrB, $bpB-$buffer, $bpB+$buffer, $bamFile, $noDuplicates);
    foreach ( @{$readRef} ) { $uniqueEntries{$_} = 1; }
    my @reads = keys %uniqueEntries;
    $readCount = scalar(@reads);

    # Count number of above reads that have unique sequence
    $readRef = &uniqueSamSequenceReads(\@reads);
    $uniqueReadCount = scalar(@{$readRef});

    # Align reads to assembly contigs using cross_match
    # Need to convert the reads to fasta format and write to a file
    my $tempFile = "/tmp/tempReads.fasta.$chrA.$bpA.$chrB.$bpB".rand();
    my $writeQualityFile = 0;

    # This returns the number of sequences written to file. It there were none
    # then can't do anything
    if ( &convertSamToFasta(\@reads, $tempFile, $writeQualityFile) == 0 ) {
        return ($readCount, $uniqueReadCount, undef, undef);
    }

    # Align reads to assembly contigs
    $crossMatchResults = &runCrossMatch($tempFile, $contigSequenceFile, $crossMatchParameters);
    $hitsToAssemblyRef = &createHitObjects($crossMatchResults);

    # Align reads to normal
    $crossMatchResults = &runCrossMatch($tempFile, $referenceSequenceFile, $crossMatchParameters);
    $hitsToNormalRef = &createHitObjects($crossMatchResults);

    # Remove the files used/created by cross_match
    unlink $tempFile;
    unlink "$tempFile.log";
    if ( -e "$tempFile.qual" ) { unlink "$tempFile.qual"; }

    return ($readCount, $uniqueReadCount, $hitsToAssemblyRef, $hitsToNormalRef);

}

sub svSpecificHits {
    # INPUT: $fastaHeader, $hitsToAssemblyRef, $hitsToNormalRef, $extend 
    #        where $extend: amount the read has to extend beyond SV breakpoint to be called spanning breakpoint
    # RETURN: Hits that align uniquely to SV assembly contig: ref to hash with key = read name; value = Hits object  
    #         Hits from above that cross SV breakpoint(s): ref to hash with key = read name; value = Hits object 
    # return ( \%specificToSv, \%crossesBreakpoints )

    my ( $fastaHeader, $hitsToAssemblyRef, $hitsToNormalRef, $extend, $maxUnalignedEndBases, $maxPercentSubs, $maxPercentIndels, $minScore) = @_;

    my ( %specificToSv, %crossesBreakpoints, $contigStart, $contigStop );

    # The breakpoints on the SV contig are encoded in the fasta header
    ($contigStart, $contigStop) = &svBreakpoints($fastaHeader);

    # See how many reads uniquely hit assembly contig and how many also cross breakpoint
    foreach my $hitName (keys %{$hitsToAssemblyRef}) {

        # Skip if this read is not aligned to the assembly contig of interest
        # The subject should equal the fasta header
        if ( $$hitsToAssemblyRef{$hitName}->subjectName() ne $fastaHeader ) { next; }

        # Skip if read does not pass the filter
        if ( !&crossMatchHitPassesFilter($$hitsToAssemblyRef{$hitName}, $maxUnalignedEndBases, $maxPercentSubs, $maxPercentIndels, $minScore ) ) { next; }

        # Skip if this read hits reference and the alignment passes filter. 
        if ( defined $$hitsToNormalRef{$hitName} && &crossMatchHitPassesFilter($$hitsToNormalRef{$hitName}, $maxUnalignedEndBases, $maxPercentSubs, $maxPercentIndels, $minScore ) ) {
            next; 
        }

        # This read uniquely hits assembly contig; it does not hit reference
        $specificToSv{$hitName} = $$hitsToAssemblyRef{$hitName};

        # Now see if hit crosses breakpoints. 
        # If the event is an insertion, only require reads to cross one of the breakpoints
        # All other events, read has to cross both breakpoints
        # The file format should probably be improved so it gives a range for both breakpoints
        if ( $fastaHeader =~ /\.INS\./ ) {
            if ( &crossMatchHitCrossesBreakpoints($$hitsToAssemblyRef{$hitName}, $contigStart-$extend, $contigStart+$extend) ||
                 &crossMatchHitCrossesBreakpoints($$hitsToAssemblyRef{$hitName}, $contigStop-$extend, $contigStop+$extend)
                ) { $crossesBreakpoints{$hitName} = $$hitsToAssemblyRef{$hitName}; }

        } else  {
            # Event is not an insertion; reads have to cross both breakpoints
            # For inversions, only one of two breakpoints is reported so the two values represents the range for the given breakpoint
            if ( &crossMatchHitCrossesBreakpoints($$hitsToAssemblyRef{$hitName}, $contigStart-$extend, $contigStop+$extend) ) {
                $crossesBreakpoints{$hitName} = $$hitsToAssemblyRef{$hitName};
            }
        }
    } # matches 'foreach my $hitName'.  Finished all hits to assembly contig '$id'

    return ( \%specificToSv, \%crossesBreakpoints);

} # end of sub


sub normalizedReadCount {
    # Input: Total number of reads from normal
    #        Total number of reads from tumor
    #        Number of SV-supporting reads from normal
    #        <optional> Percent of normal sample that consists of tumor reads
    # Return: Normalized number of SV-supporting reads from normal sample after
    #         correcting for sample purity and differences in coverage

    my ($totalNormalReads, $totalTumorReads, $normalSvReads, $percentTumorInNormal) = @_;
    if ( $totalNormalReads == 0 || $normalSvReads == 0 ) { return 0; }

    if ( !defined $percentTumorInNormal ) { $percentTumorInNormal = 0; }

    my $adjustedNormalSvReads = $normalSvReads - ($percentTumorInNormal * $normalSvReads)/100;
    my $normalizedNormalSvReads = $adjustedNormalSvReads * $totalTumorReads/$totalNormalReads;

    # Round to nearest integer
    $normalizedNormalSvReads = int($normalizedNormalSvReads + 0.5); 

    return $normalizedNormalSvReads;
}


sub displayAlignments {
    # For debugging
    # Input: ref to hash with key = read name; value = Hits object 
    #        (There can be two inputs; one for read to SV contig, one for read to reference sequence
    # Output: print to STDOUT the alignment lines for each read

    my ( $hits1Ref, $hits2Ref ) = @_;
    if ( !defined $hits2Ref ) { $hits2Ref = (); }

    foreach my $name (keys %{$hits1Ref}) {
        print $$hits1Ref{$name}->alignmentLine(), "\n";
        if ( defined $$hits2Ref{$name} ) {
            print $$hits2Ref{$name}->alignmentLine(), "\n";
        }
    }
}




sub normalAlleleCount {
    # return ($leftCrossingBreakpoint, $leftTotalReads, $rightCrossingBreakpoint, $rightTotalReads);
    # need to decide how to deal with ambiguity....
    # right now, assuming that the outer boundaries of event are chosen


    my ($fastaHeader, $hitsToAssemblyRef, $bamFile) = @_;

    my ($contigStart, $contigStop, $ambiguity, $chrA, $bpA, $chrB, $bpB, $event, $subjectStart,
        $subjectStop, $noDuplicates, $readRef, $read, $leftTotalReads, $leftCrossingBreakpoint, 
        $rightTotalReads, $rightCrossingBreakpoint, $name, $flag, $cigar, );

    $noDuplicates = 1;
    my $maxMismatch = 3;
    my $maxSoftmask = 3;

    $rightTotalReads = $rightCrossingBreakpoint = $leftTotalReads = $leftCrossingBreakpoint = 0;

    # Get SV breakpoints in reference from fasta header.
    if ( $fastaHeader =~ /Var\:(\w+)\.(\d+)\.(\w+)\.(\d+)\.(\w+)\./ ) {
        $chrA = $1; $bpA = $2; $chrB = $3; $bpB = $4; $event = $5;
    } else {
        confess "Unexpected format for fasta header '$fastaHeader'";
    }

    # Get breakpoint(s) in SV contig encoded in the fasta header.  The ambiguity in the SV breakpoint position is 
    # the difference between $contigStart and $contigStop.
    ($contigStart, $contigStop) = svBreakpoints($fastaHeader);
    $ambiguity = $contigStop - $contigStart - 1;


    # Get non-duplicated reads around breakpoints
    # Then see which of them span the reference breakpoint regions

    $readRef = &getReads($chrA, $bpA, $bpA, $bamFile, $noDuplicates);
    $leftTotalReads = scalar( @{$readRef} );
    foreach $read ( @{$readRef}  ) {
        chomp $read;

        # Get name that was used for this read as a fasta sequence in cross_match
        ($name, $flag, undef, undef, undef, $cigar) = split /\s+/, $read;
        $name = "$name.$cigar.$flag";

        # Don't count this read if it has an alignment to SV contig that passes filters
        if ( defined $$hitsToAssemblyRef{$name} && &crossMatchHitPassesFilter($$hitsToAssemblyRef{$name}) ) {
            next;
        }

        # If it passes the filters for BWA-aligned reads and it crosses breakpoint + ambiguity, count it as supporting reference allele
        if ( samReadPassesFilter($read, $maxMismatch, $maxSoftmask) && 
             &samReadCrossesBreakpoints($read, $bpA, $bpA+$ambiguity) ) { 
            $leftCrossingBreakpoint++; 
        }
    }

    $readRef = &getReads($chrB, $bpB, $bpB, $bamFile, $noDuplicates);
    $rightTotalReads = scalar( @{$readRef} );
    foreach $read ( @{$readRef}  ) {
        chomp $read;

        # Get name that was used for this read as a fasta sequence in cross_match
        ($name, $flag, undef, undef, undef, $cigar) = split /\s+/, $read;
        $name = "$name.$cigar.$flag";

        # Don't count this read if it has an alignment to SV contig that passes filters
        if ( defined $$hitsToAssemblyRef{$name} && &crossMatchHitPassesFilter($$hitsToAssemblyRef{$name}) ) {
            next;
        }

        # If it passes the filters for BWA-aligned reads and it crosses breakpoint - ambiguity, count it as supporting reference allele
        if ( samReadPassesFilter($read, $maxMismatch, $maxSoftmask) && 
             &samReadCrossesBreakpoints($read, $bpB-$ambiguity, $bpB) ) { 
            $rightCrossingBreakpoint++; 
        }
    }

    return ($leftCrossingBreakpoint, $leftTotalReads, $rightCrossingBreakpoint, $rightTotalReads);

}






sub sequenceSpecificHits {
    # generic version of svSpecificHits that can be used to get reads specific for reference or SV contig
    # INPUT: $hitsToContigRef, $hitsToCompetingContigRef, $contigStart, $contigStop, $extend 
    # $contigStart, $contigStop are positions on contig that must be spanned.  The coordinates are relative
    # to the contig (so if the contig is genomic sequence, need to convert genomic coordinates to what is in the fasta
    # A sequence-specific read is not allowed to align to a competing contig

}




sub alleleCount  {
    # testing

    # Hacked up version that will look at multiple alignments of a given read to the reference.
    # This is because the two breakpoint fasta sequences could overlap for breakpoints that are close

    my ($fastaHeader, $hitsToAssemblyRef, $hitsToNormalRef, $extend, $crossMatchResults) = @_;

    my ($contigStart, $contigStop, $ambiguity, $chrA, $bpA, $chrB, $bpB, $event, $svBreakPointReadsRef,
        $subjectName, $subjectChr, $subjectStart, $subjectStop, %crossesRight, %crossesLeft, );

    # Get strict list of reads that cross SV breakpoint
    ## should make  &svSpecificHits generic so can use for both reference hits and
    #  SV hits
    (undef, $svBreakPointReadsRef) = &svSpecificHits($fastaHeader, $hitsToAssemblyRef, $hitsToNormalRef, $extend);

    # Get breakpoint(s) in SV contig encoded in the fasta header.  The ambiguity in the SV breakpoint position is 
    # the difference between $contigStart and $contigStop.
    ($contigStart, $contigStop) = svBreakpoints($fastaHeader);
    $ambiguity = $contigStop - $contigStart - 1;

    # Get SV breakpoints in reference from fasta header.
    if ( $fastaHeader =~ /Var\:(\w+)\.(\d+)\.(\w+)\.(\d+)\.(\w+)\./ ) {
        $chrA = $1; $bpA = $2; $chrB = $3; $bpB = $4; $event = $5;
    } else {
        confess "Unexpected format for fasta header '$fastaHeader'";
    }

    foreach my $line ( @{$crossMatchResults} ) {
        chomp $line;

        if ( $line !~ /$AlignmentLine/ ) { next; }


        my $hit = new Hits;
        $hit->addCrossMatchLine($line);
        my $query = $hit->queryName();

        # Skip if read does not pass the filter
        # Skip if this read hits SV contig and the alignment to SV contig passes filter. 
        if ( !&crossMatchHitPassesFilter($hit) ||
             ( defined $$hitsToAssemblyRef{$query} && &crossMatchHitPassesFilter($$hitsToAssemblyRef{$query}) )
            ) {
            next; 
        }
        # See if the read crosses one of the breakpoints.  Skip if the breakpoint is not within subject
        # boundaries encoded in the subject name (fasta header of reference sequence)
        $subjectName = $hit->subjectName();
        if ( $subjectName =~ /(\w+)\.(\d+)\.(\d+)/ ) {
            $subjectChr = $1; $subjectStart = $2; $subjectStop = $3;
        } else {
            confess "Unexpected format for reference contig fasta header (i.e. subject of Hits object): '$subjectName'";
        }
        # Look at left breakpoint
        if ( $subjectChr eq $chrA && $bpA > $subjectStart && $bpA < $subjectStop ) {
            # Range the read must cover is $bpA, $bpA+$ambiguity
            # Need to convert genomic coordinates to coordinates used in alignment to reference
            # (subtract $subjectStart and add 1)
            if ( &crossMatchHitCrossesBreakpoints($hit, $bpA-$subjectStart+1, $bpA+$ambiguity-$subjectStart+1) ) {
                $crossesLeft{$query} = $hit;
            }

        }

        # Look at right breakpoint
        if ( $subjectChr eq $chrB && $bpB > $subjectStart && $bpB < $subjectStop ) {
            # Range the read must cover is $bpB-$ambiguity, $bpB
            # Need to convert genomic coordinates to coordinates used in alignment to reference
            # (subtract $subjectStart and add 1)
            if ( &crossMatchHitCrossesBreakpoints($hit, $bpB-$ambiguity-$subjectStart+1, $bpB-$subjectStart+1) ) {
                $crossesRight{$query} = $hit;
            }
        }
    } # matches 'foreach $line'

    return($svBreakPointReadsRef, \%crossesLeft, \%crossesRight);

} # end of sub




sub BAK_alleleCount {
    # INPUT: $fastaHeader, $hitsToAssemblyRef, $hitsToNormalRef, $extend 
    #        where $extend: amount the read has to extend beyond SV breakpoint to be called spanning breakpoint
    # RETURN: Hits that support SV contig allele: ref to hash with key = read name; value = Hits object  
    #         Hits that support left Reference breakpoint: ref to hash with key = read name; value = Hits object  
    #         Hits that support right Reference breakpoint: ref to hash with key = read name; value = Hits object
    # return ( \%svReads, \%leftBreakpoint, \%rightBreakpoint )
    # 
    # Assume the alignments are to reference made by &getBuild36ReferenceSequences 
    # In the Hits objects, the subject name will be '$chr.$start.$stop' so can convert to genomic
    # coordinates from subject coordinates by adding '$start'  (and add 1 or subtract 1 ???)

    # For assembly contig breakpoint ambiguity, require reads cross entire region of ambiguity
    # For reference, look at left and right breakpoints separately and take a range of single 
    # breakpoint based on ambiguity.  
    # How the ambiguity in SV contig breakpoint is handled depends on the type of event and how
    # Ken reports genomic breakpoints.  For CTX, it will also depend on the orientation of alignments
    # So could add or subtract the ambiguity to the breakpoint reported in the fasta header
    # If Ken reports the 'outer' breakpoints, then the range for the left will be $bpA, $bpA + $ambiguity
    # and range for right is $bpB - $amgiguity, $bpB
    #
    # For now, just add ambiguity to left and subtract from right of reference


    #  OOPS: the method for making Hits objects only returns one hit for a given read.  
    # When taking +/- 500 bp around each breakpoint, there will be fasta entries that overlap.
    # But will not get both hits to overlapping sequence.  Will arbitrarily choose one of the hits.
    #     That is not a problem.  Just parse the subject header to see if the breakpoint
    #     falls within it.  It doesn't matter which of the two fasta entries was chosen


    my ($fastaHeader, $hitsToAssemblyRef, $hitsToNormalRef, $extend) = @_;

    my ($contigStart, $contigStop, $ambiguity, $chrA, $bpA, $chrB, $bpB, $event, $svBreakPointReadsRef,
        $subjectName, $subjectChr, $subjectStart, $subjectStop, %crossesRight, %crossesLeft, );

    # Get strict list of reads that cross SV breakpoint
    ## should make  &svSpecificHits generic so can use for both reference hits and
    #  SV hits
    (undef, $svBreakPointReadsRef) = &svSpecificHits($fastaHeader, $hitsToAssemblyRef, $hitsToNormalRef, $extend);


    # Get breakpoint(s) in SV contig encoded in the fasta header.  The ambiguity in the SV breakpoint position is 
    # the difference between $contigStart and $contigStop.
    ($contigStart, $contigStop) = svBreakpoints($fastaHeader);
    $ambiguity = $contigStop - $contigStart - 1;

    # Get SV breakpoints in reference from fasta header.
    if ( $fastaHeader =~ /Var\:(\w+)\.(\d+)\.(\w+)\.(\d+)\.(\w+)\./ ) {
        $chrA = $1; $bpA = $2; $chrB = $3; $bpB = $4; $event = $5;
    } else {
        confess "Unexpected format for fasta header '$fastaHeader'";
    }


    # See how many reads uniquely hit reference and also cross breakpoint
    foreach my $hitName (keys %{$hitsToNormalRef}) {

        # Skip if read does not pass the filter
        # Skip if this read hits SV contig and the alignment to SV contig passes filter. 
        if ( !&crossMatchHitPassesFilter($$hitsToNormalRef{$hitName}) ||
             ( defined $$hitsToAssemblyRef{$hitName} && &crossMatchHitPassesFilter($$hitsToAssemblyRef{$hitName}) )
            ) {
            next; 
        }

        # See if the read crosses one of the breakpoints.  Skip if the breakpoint is not within subject
        # boundaries encoded in the subject name (fasta header of reference sequence)
        $subjectName = $$hitsToNormalRef{$hitName}->subjectName();
        if ( $subjectName =~ /(\w+)\.(\d+)\.(\d+)/ ) {
            $subjectChr = $1; $subjectStart = $2; $subjectStop = $3;
        } else {
            confess "Unexpected format for reference contig fasta header (i.e. subject of Hits object): '$subjectName'";
        }

        # Look at left breakpoint
        if ( $subjectChr eq $chrA && $bpA > $subjectStart && $bpA < $subjectStop ) {
            # Range the read must cover is $bpA, $bpA+$ambiguity
            # Need to convert genomic coordinates to coordinates used in alignment to reference
            # (subtract $subjectStart and add 1)
            if ( &crossMatchHitCrossesBreakpoints($$hitsToNormalRef{$hitName}, $bpA-$subjectStart+1, $bpA+$ambiguity-$subjectStart+1) ) {
                $crossesLeft{$hitName} = $$hitsToNormalRef{$hitName};
            }

        }

        # Look at right breakpoint
        if ( $subjectChr eq $chrB && $bpB > $subjectStart && $bpB < $subjectStop ) {
            # Range the read must cover is $bpB-$ambiguity, $bpB
            # Need to convert genomic coordinates to coordinates used in alignment to reference
            # (subtract $subjectStart and add 1)
            if ( &crossMatchHitCrossesBreakpoints($$hitsToNormalRef{$hitName}, $bpB-$ambiguity-$subjectStart+1, $bpB-$subjectStart+1) ) {
                $crossesRight{$hitName} = $$hitsToNormalRef{$hitName};
            }
        }

    } # matches 'foreach my $hitName'   

    return($svBreakPointReadsRef, \%crossesLeft, \%crossesRight);

} # end of sub


sub fractionComplemented {
    # INPUT: ref to hash of Hits.pm objects
    # RETURN: fraction of these hits where the subject is complemented
    #   When converting SAM format to fasta, the reads are not reverse complemented
    #   depending on flag.  Get flag from name to see if it should be complemented

    my $hitsRef = $_[0];
    my ( $total, $complemented, $fraction, $flag, $readToReferenceReversed, );
    foreach my $hit ( keys %{$hitsRef} ) {
        $readToReferenceReversed = 0;
        # Get flag that is appended to read name in sub convertSamToFasta
        if ( $$hitsRef{$hit}->queryName() =~ /\.(\d+)$/ ) { $flag = $1; } else { confess "Did not get flag for Hit object '", $$hitsRef{$hit}->queryName(), "'"; }
        if ( $flag & 0x0010 ) { $readToReferenceReversed = 1; }
        if ( $$hitsRef{$hit}->isComplemented() != $readToReferenceReversed ) { $complemented++; }
        $total++;
    }

    if ( $total == 0 ) { return 0; }
    $fraction = $complemented/$total;

    # Round off fraction to nearest 0.001
    $fraction += 0.0005;
    if ( $fraction =~ /(\d+\.\d{3})/ ) { $fraction = $1; }

    return $fraction;
}

sub svBreakpoints {
    # INPUT: fasta header
    # RETURN: SV breakpoints start and stop

    my $fastaHeader = $_[0];
    my ($contigStart, $contigStop);

    # The breakpoints on the SV contig are encoded in the fasta header
    if ( $fastaHeader =~ /Ins\:(\d+)\-(\-|\d+)/ ) {
        $contigStart = $1; $contigStop = $2;
    } else {
        confess "Unexpected format for fasta header.  Did not get breakpoints on SV contig: '$fastaHeader'";
    }
    ( defined $contigStart && defined $contigStop ) ||
        confess "Did not get breakpoints on SV contig from '$fastaHeader'";
    if ( $contigStop eq "-" ) { $contigStop = $contigStart; }

    return ($contigStart, $contigStop);
}


1;
