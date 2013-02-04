package Genome::Model::Tools::Sv::SvFusionGenes;

# fusionGenes.pl
#
# Requires input file to be in BreakDancer format
#
#   Require both breakpoints to be in introns of different genes (for now)
#   This should work for all events as long as we know relative orientations of the junction sequences
# 
# General Tactics:
# Get BreakDancer line, get first and second breakpoint
# Get every transcript that crosses breakpoint
# Fusion is possible if there are different genes for transcripts
# Require that both breakpoints are in introns
# Need relative oriention of 'fused' sequences
#     +- Both increasing (this will always be case for DEL)
#     -+ Both decreasing
#     ++ Increasing A, decreasing B
#     -- Decreasing A, increasing B
# Based on above, look at strands the genes are on and determine if a fusion transcript is possible
# Look for fusion transcripts without frameshift



use strict;
use warnings;
use Carp;
use Genome;

class Genome::Model::Tools::Sv::SvFusionGenes {
    is => 'Genome::Model::Tools::Sv',
    has => [
        input_file => {
            type => 'String',
            doc => 'Input files in breakdancer format',
            is_optional => 0,
        },
        output_file => {
            type => 'String',
            doc => 'output sv annotation file',
        },
        annotation_build_id => {
            is => 'Text',
            doc => 'Annotation build you want to use',
            #default => '102549985',
            default => '131184146',
        },
        annotation_build => {
            is => "Genome::Model::Build::ImportedAnnotation",
            id_by => 'annotation_build_id',
        },
    ],
};


sub execute {
    my $self  = shift;
    my $build = $self->annotation_build;
    unless ($build){
        die "build not defined";
    }

    $self->processFile;
}


=head

NEED TO USE THE FASTA HEADER
always before ,Ins:
14.13,tumor14,Var:14.67672785.14.67673511.DEL.700.+-,Ins:560-586,Length:1341,KmerCoverage:29.01,Strand:+,Assembly_Score:1257.67,PercNonRefKmerUtil:0,TIGRA

CTX ++
CTX ++,++
CTX ++,--
CTX +-
CTX +-,+-
CTX +-,-+
CTX +-,--
CTX -+
CTX -+,+-
CTX -+,-+
CTX -+,--
CTX --
CTX --,++
CTX --,--
DEL +-
DEL +-,+-
INS +-
INS +-,+-
INV ++
INV --
INV --,--
ITX +-
ITX +-,+-
=cut



sub processFile {
    my ( @entireFile, $line, $bdRef, $chrA, $bpA, $chrB, $bpB, $event, $patient, $aTranscriptRef, $bTranscriptRef, $transcriptRef,
	 $eventId, $geneA, $transcriptA, $orientationA, $subStructureA, $geneB, $transcriptB, $orientationB, $subStructureB,
	 $deletedGeneHashRef, $deletedGenes, $inCommonRef, $transcript, $junctionOrientation,
	);

    my $in_fh  = Genome::Sys->open_file_for_reading($self->input_file) or die "Failed to open input file\n";
    my $out_fh = Genome::Sys->open_file_for_writing($self->output_file) or die "Failed to open output file\n";
    
    while (my $line = $in_fh->getline) { 
	    chomp $line;
	    next if $line =~ /^\s*$/ or $line =~ /\#/;
	
	    ($eventId, $chrA, $bpA, undef, $chrB, $bpB, undef, $event, $junctionOrientation ) = split /\s+/, $line;
	    $geneA = $transcriptA = $orientationA = $subStructureA = $geneB = $transcriptB = $orientationB = $subStructureB = $deletedGenes = "N/A";
	    $inCommonRef = $deletedGeneHashRef = undef;

	    # get junction orientation
	    if ( $line =~ /([\+\-]*)\,Ins\:/ ) { 
            $junctionOrientation = $1; 
        }

	    my @cols = split /\t/, $line;
       	
	    # get all transcripts crossing each breakpoint
	    $aTranscriptRef = allTranscripts($chrA, $bpA, $build);
	    $bTranscriptRef = allTranscripts($chrB, $bpB, $build);
	
	    # If there are no transcripts crossing either breakpoint, then you are done
	    next if !defined $aTranscriptRef || scalar @{$aTranscriptRef} == 0 || !defined $bTranscriptRef || scalar @{$bTranscriptRef} == 0;
	    
	    # See if there are different genes crossing the breakpoints
	    if ( fusionGenes($aTranscriptRef, $bTranscriptRef) ) {
	        my $aGeneNameRef = geneNames($aTranscriptRef);
	        my $bGeneNameRef = geneNames($bTranscriptRef);
	        my $fusionProteinRef = getAllInFrameFusions($aTranscriptRef, $bpA, $bTranscriptRef, $bpB, $junctionOrientation);
	        if ( defined $fusionProteinRef && scalar(keys%{$fusionProteinRef}) >= 1 ) {
		        $out_fh->print("\n$line\t[");
                map{$out_fh->print("$_ ")}sort keys %{$aGeneNameRef};
		        $out_fh->print("]\t[");
		        map{$out_fh->print("$_ ")}sort keys %{$bGeneNameRef};
                $out_fh->print("] \n");

		        for my $seq (keys %$fusionProteinRef) { 
		            for my $mRNA ( keys %{$fusionProteinRef->{$seq}} ) {
			            for my $mRNA_withCoordinates ( keys %{$fusionProteinRef->{$seq}->{$mRNA}} ) {
			                $out_fh->print("\t$seq\n\t$mRNA\n\t$mRNA_withCoordinates\n"); 
                        }
                    }
                }
            }
 
        }
    }

    $out_fh->close;
}


sub getAllInFrameFusions {
    # Do all pairwise combinations of transcripts
    # Require both breakpoints to be in an intron (for now)
    # Check to see if relative orientation of sequences and genes make a fusion on one strand
    # Only include the fusions that are still in frame after putting the exons together

    my ($aTranscriptRef, $bpA, $bTranscriptRef, $bpB, $junctionOrientation) = @_;
    my ( %allFusions, $mRNA_a, $mRNA_b, $substructureA,  $substructureB, $firstSecond, $fusionMessage, $seqObj, 
	 $protObj, $protSeq,  $fivePrimeMessage, $threePrimeMessage, $rnaWithCoordinates_a, $rnaWithCoordinates_b,
	 $fivePrimeMessageWithCoordinates, $threePrimeMessageWithCoordinates, 
	);

    for my $aTranscript ( @{$aTranscriptRef} ) {
	    # Confirm breakpoint is in an intron
	    $substructureA = substructureWithBreakpoint($aTranscript, $bpA);
	    next if $substructureA !~ /intron/;

	    for my $bTranscript ( @{$bTranscriptRef} ) {
	        # Confirm breakpoint is in an intron
	        $substructureB = substructureWithBreakpoint($bTranscript, $bpB);
	        next if $substructureB !~ /intron/;

	        # Determine which transcript is 5' and which is 3' of fused message based on
	        # relative orientation of sequences and strand that each gene is on
	        $firstSecond = firstSecond($junctionOrientation, $aTranscript, $bTranscript);
	        # If transcripts are not going in same direction, there is no in-frame fusion possible and '$firstSecond' is returned undef
	        next unless defined $firstSecond;

	        # If both substructures are exons, get assembly contig and write to a file
	        # Need to check to see if event is CTX or other event and pass correct file to be read
	        # (i.e. need both assembly files as input so this won't work if a lot of different patient samples
	        #  are in file.  Will need to separate out patients or do file of files for assembly *fasta files.
	        #  Will also need to code an identifier for each patient.  Gets too complicated.


	        if ( $firstSecond eq "AB" ) {
		        ($mRNA_a, $rnaWithCoordinates_a) = processedMessage($aTranscript, $bpA, "start");
		        ($mRNA_b, $rnaWithCoordinates_b) = processedMessage($bTranscript, $bpB, "end");
		        $rnaWithCoordinates_a = "<".$aTranscript->gene_name."_". $aTranscript->transcript_name."> ".$rnaWithCoordinates_a;
		        $rnaWithCoordinates_b = "<".$bTranscript->gene_name."_". $bTranscript->transcript_name."> ".$rnaWithCoordinates_b;
            } 
            elsif ( $firstSecond eq "BA" ) {
		        ($mRNA_a, $rnaWithCoordinates_a) = processedMessage($bTranscript, $bpB, "start");
		        ($mRNA_b, $rnaWithCoordinates_b) = processedMessage($aTranscript, $bpA, "end");
		        $rnaWithCoordinates_a = "<".$bTranscript->gene_name."_". $bTranscript->transcript_name."> ".$rnaWithCoordinates_a;
		        $rnaWithCoordinates_b = "<".$aTranscript->gene_name."_". $aTranscript->transcript_name."> ".$rnaWithCoordinates_b
            }
	        else {
		        confess "Unexpected \$firstSecond: '$firstSecond'";
            }

	        $fusionMessage = $mRNA_a.$mRNA_b;
	        next if $fusionMessage eq "";
	        $seqObj  = new_sequence($fusionMessage);
	        $protObj = $seqObj->translate;
	        $protSeq = $protObj->seq;

	        # See that the fusion protein is in frame
	        if ( $protSeq !~ /\w\*\w/ && $protSeq =~ /\w\*$/ ) {
		        my $messageWithCoordinates = $rnaWithCoordinates_a. " | ".$rnaWithCoordinates_b;
		        my $message = $messageWithCoordinates;

		        # This gets rid of all of the extra annotation of the sequence
		        $message =~ s/\[\d+\]\s+\[\d+\]//g; 
                $message =~ s/\[\d+\]//g; 
                $message =~ s/\s*//g; 
                $message =~ s/\<\w+\_\w+\>//g;
		        $message =~ s/\|/ | /; 
		        $allFusions{$protSeq}{$message}{$messageWithCoordinates} = 1;  # put both protein sequences here, separated by '|' or ....
		
		        #  need to get translation of each end, but a codon could be split... first half could have different frame than second half
		        my ($halfSeq, $halfProtObj, $halfProtSeq);
		        $halfProtSeq = "--";
		        if ( $mRNA_a ne "" ) {
		            $halfSeq = new_sequence($mRNA_a); $halfProtObj = $halfSeq->translate; $halfProtSeq = $halfProtObj->seq;
                }
		        $out_fh->print("\n\n\$mRNA_a \t $halfProtSeq \n");

		        $halfProtSeq = "--";
		        if ( $mRNA_b ne "" ) {
		            $halfProtSeq = inFrameSequence($mRNA_b);
		        }
		        $out_fh->print("\$mRNA_b \t $halfProtSeq \n\n");
            }
        }
    }
    return \%allFusions;
}

sub inFrameSequence {
    my $message = shift;
    my ( $seqObj, $protObj, $protSeq  ); 
    
    for my $i (1..3) {
	    $seqObj = new_sequence($message);
	    $protObj = $seqObj->translate;
	    $protSeq = $protObj->seq;
	    if ( $protSeq !~ /\w\*\w/ ) { 
            return $protSeq; 
        }
	    # Change reading frame
	    if ( $message =~ /\w(\w+)/ ) { 
            $message = $1; 
        }
    }
    return "None in frame";
    
}

sub processedMessage {
    # Return message up to given position
    # Need to specify if want start or end of gene....
    # For now, assume the position is within an intron
    

    # For 5' part of message, do not include utr_exon in translated message
    # For 3' part of message, need to include utr_exon that occur after breakpoint and before first cds_exon

    my ( $transcript, $position, $startOrEnd ) = @_;
    my $gene = $transcript->gene_name; my $transcriptName = $transcript->transcript_name;
    my @ss   = $transcript->ordered_sub_structures;

    my $mRNA = "";
    my $rnaWithCoordinates = "";
    my $reading = 0;
    my $includeUtrExons = 1;  # for 3'
    if ( $startOrEnd eq "start" ) { 
	    $includeUtrExons = 0;  # for 5'
	    $reading = 1; 
    }
    my ($left, $right);
    for my $sub_structure ( @ss ) {
	    # See if the position is between the substructure start and stop
	    if ( $position >= $sub_structure->structure_start && $position <= $sub_structure->structure_stop ) { 
	    # If this was the 3' part of message (i.e. 'start' of message), then quit reading substructures after the breakpoint
	        if ( $startOrEnd eq "start" ) { 
		        $reading = 0; 
            }

	        # If this is the 5' part of message, then start reading substructures after breakpoint
	        if ( $startOrEnd eq "end" ) { 
		        $reading = 1; 
            }
        }
	
	    # Just keep going if the substructure is not part of fused message
	    next unless $reading;
 
	    # 'flank' does not have sequence in database and don't want 'intron'
	    next if $sub_structure->structure_type eq "intron" or $sub_structure->structure_type eq 'flank';
	
    	# The message to be translated has cds_exon and if part of the 5' has utr_exon that are beginning of that transcript
	    if ( $sub_structure->structure_type eq "cds_exon" ) {	     
	        $mRNA .= $sub_structure->nucleotide_seq; 
	        # Only include utr_exons until the first cds_exon
	        $includeUtrExons = 0; 
	    } 
        elsif ( $sub_structure->structure_type eq "utr_exon" and $includeUtrExons ) {	     
	        $mRNA .= $sub_structure->nucleotide_seq; 
        }


	    # This is the full message with coordinates.  It is used to design primers
	    $left  = " <".$sub_structure->structure_type. $sub_structure->ordinal."> [".$sub_structure->structure_start."]";
	    $right = "[".$sub_structure->structure_stop."]";
	    if ( $transcript->strand == 1 ) {
	        $rnaWithCoordinates = $rnaWithCoordinates.$left.$sub_structure->nucleotide_seq.$right;
	    } 
        elsif ($transcript->strand == -1 ) {
	        $rnaWithCoordinates = $rnaWithCoordinates.$right.$sub_structure->nucleotide_seq.$left;
	    } 
        else {
	        confess "Unexpected strand value";
        }		
    }
    return ($mRNA, $rnaWithCoordinates);
}


sub firstSecond {
    # Returns either "AB" or "BA" indicating which gene is 5' and which is 3'
    my ( $junctionOrientation, $aTranscript, $bTranscript ) = @_;

    my ( $aStrand, $bStrand, $firstSecond, );
    # Following returns -1 or 1
    $aStrand = $aTranscript->strand;
    $bStrand = $bTranscript->strand;

    if ( $aStrand != 1 && $aStrand != -1 ) {
        die "\nDid not get strand for aTranscript: $aTranscript"; 
    }
    if ( $bStrand != 1 && $bStrand != -1 ) { 
        die "\nDid not get strand for bTranscript: $bTranscript"; 
    }
    
    if ( $junctionOrientation eq "+-" ) {
	    # Both sequences increasing.  Transcripts must be on same strand to have in frame fusion.
	    #  if $aStrand and $bStrand have opposite strands, no in frame fusion
	    if ( $aStrand == 1 && $bStrand == 1 ) { 
            $firstSecond = "AB"; 
        }
	    if ( $aStrand == -1 && $bStrand == -1 ) { 
            $firstSecond = "BA"; 
        }
    }
    elsif ( $junctionOrientation eq "-+" ) {
	    # Both sequences decreasing.  Transcripts must be on same strand to have in frame fusion.
	    #  if $aStrand and $bStrand have opposite strands, no in frame fusion
	    if ( $aStrand == -1 && $bStrand == -1 ) { 
            $firstSecond = "AB"; 
        }
	    if ( $aStrand == 1 && $bStrand == 1 ) { 
            $firstSecond = "BA"; 
        }
    }   
    elsif ( $junctionOrientation eq "++" ) {
	    # Sequence A increasing, sequence B decreasing. Transcripts must be on opposite strands to have in frame fusion
	    if ( $aStrand == 1 && $bStrand == -1 ) { 
            $firstSecond = "AB"; 
        }
	    if ( $aStrand == -1 && $bStrand == 1 ) { 
            $firstSecond = "BA"; 
        }
    } 
    elsif ( $junctionOrientation eq "--" ) {
	    # Sequence A decreasing, sequence B increasing. Transcripts must be on opposite strands to have in frame fusion
	    if ( $aStrand == -1 && $bStrand == 1 ) { 
            $firstSecond = "AB"; 
        }
	    if ( $aStrand == 1 && $bStrand == -1 ) { 
            $firstSecond = "BA"; 
        }
    } 
    else {
	    confess "Unexpected junction orientation: '$junctionOrientation'";
    }

    return $firstSecond;
}


sub allTranscripts {
    # Return all transcripts spanning the given position
    my ( $chr, $position, $build ) = @_;
    
    my ( $transcriptIterator, $transcript, @transcripts, );
    $transcriptIterator = $build->transcript_iterator(chrom_name => $chr); 
    die "transcriptIterator not defined" unless defined $transcriptIterator;
    while ( $transcript = $transcriptIterator->next ) {
	    if ( $position >= $transcript->transcript_start and $position <= $transcript->transcript_stop ) {
	        push @transcripts, $transcript;
	    }
    }
    
    return \@transcripts;
}


sub fusionGenes {
    # Input: ref to transcripts crossing each breakpoint
    # Return: 1 if there are different genes crossing the breakpoints
    my ( $aTranscriptRef, $bTranscriptRef ) = @_;
    my $aGeneNameRef = geneNames($aTranscriptRef);
    my $bGeneNameRef = geneNames($bTranscriptRef);

    for my $geneA (keys %$aGeneNameRef) {
	    next if $geneA eq "";
        return 1 unless defined $bGeneNameRef->{$geneA};
    }
    return 0;
}
    

sub geneNames {
    # Return ref to hash of gene names for the list of transcripts
    my $transcriptRef = shift;
    my %geneNames;
    for my $transcript ( @$transcriptRef ) {
	    my $name = $transcript->gene_name;
	    if ( defined $name ) { 
            $geneNames{$name} = 1; 
        }
    }
    return \%geneNames;
}

sub substructureWithBreakpoint {
    # Return substructure with breakpoint.  The name of substructure is the type
    # and the number of coding exons before the substructure (does not utr_exon as coding)
    # type.cds_exons_before
    
    my ($transcript, $position) = @_;
    my @subStructures = $transcript->ordered_sub_structures;

    for my $structure (@subStructures) {
	    if ( $position >= $structure->{structure_start} and $position <= $structure->{structure_stop} ) {
	        return $structure->{structure_type}.$structure->{cds_exons_before};
        }
    }
    confess "'$position' is not in transcript";
}




