package Genome::Model::Tools::Annotate::Sv::FusionTranscripts;

use strict;
use warnings;
use Carp;

use File::Basename;
use Bio::Perl;
use Genome;

class Genome::Model::Tools::Annotate::Sv::FusionTranscripts {
    is => 'Genome::Model::Tools::Annotate::Sv::Base',
    doc => "Determine whether SV may cause transcript fusions",
    has_input => [
        fusion_output_file => {
            type => 'String',
            doc  => 'output fusion annotation file',
        },
    ],
};


#The example output of MergeAssembledCallsets.pl, which is used as input file to this command
#2.21    2       89461280        89461280        2       89495566        89495566        DEL     +-      34283   34283   tumor2  1215    NA      NA      NA
#7.28    7       57808537        57808537        7       57808690        57808690        DEL     +-,+-   154     156     tumor7,normal7  1198,626        NA      NA      NA
#If there are two in tumor/normal, xxxx.svs.merge.fasta will use the one with higher asmscore; if both get the same asmscore, it just 
#uses the first one (before ,).  But how to deal with the conflict of orientation of two like (++,-+) ? This won't be a problem for
#somatic event because there is only one orientation (for tumor).

sub help_detail {
    return "This code predicts the potential transcript fusion between two breakpoints of the different 
    genes and same genes. It is based on John Wallis's fusionGenes.pl with the addition of in-frame fusion for 
    the same gene. Below is the general tactics commented in the original script:

    1. Get sv line, get first and second beakpoint
    2. Get every transcript that crosses breakpoint
    3. Fusion is possible if there are different genes for transcripts
    4. Require that both breakpoints are in introns
    5. Need relative oriention of 'fused' sequences
         +- Both increasing (this will always be case for DEL)
         -+ Both decreasing
         ++ Increasing A, decreasing B
         -- Decreasing A, increasing B
    6. Based on above, look at strands the genes are on and determine if a fusion transcript is possible
    7. Look for fusion transcripts without frameshift
    
For the case of the same gene fusion, two fused transcripts must be the same and the number difference of two exon ordinals must be > 1 (like exon3 -> exon5).

In the fusion-transcripts output, each fusion event has a set of results formatted as following mock example:

    \$mRNA_a   PPSGLLTPLHLEGLVQFQDVSFAYPNRPDVLVLQ
    \$mRNA_b   TTIMAVEF

    6    32815568    6    32816364    INV    [PSMB9 TAP1]|[PSMB9 TAP1]
    PPSGLLTPLHLEGLVQFQDVSFAYPNRPDVLVLQ \| TTIMAVEF
    <TAP1_ENST00000354258> [32821755]GAGAAGG<utr_exon1> [32821594][32821593]ATGGCTGAGC<cds_exon1> 
    [32820816][32820279]GCCAG <cds_exon2> [32820165] \| <PSMB9_ENST00000395330><cds_exon2> [32823924]
    ATGGCAGTGGAGTTTGACGGGGGCGTTGTGATGGGTTCTGATTCCCGAGTGTCTGCAGG[32823982] <utr_exon6> [32827310]
    AAACTCTCTAGGGCCAAAA[32827362]


* \$mRNA_a shows the amino acid sequence by transcript on left side of fusion point, while \$mRNA_b shows the amino acid sequence by transcript on right side of fusion point.

** The headers of '6    32815568    6    32816364    INV    [PSMB9 TAP1]|[PSMB9 TAP1]' are: chromosomeA  breakpointA  chromosomeB breakpointB  eventType  FusionGenes

*** 'PPSGLLTPLHLEGLVQFQDVSFAYPNRPDVLVLQ | TTIMAVEF' is the fused amino acid sequence and '|' is the fusion point.

**** The last part is the detailed sequence with coordinates of each substructure. It starts with <GeneName_TranscriptName> on left end , then coordinates and sequences and names of substructures (only list utr_exon and cds_exon) all the way to the fusion point '|', then <GeneName_TranscriptName> on right side of fusion point, then coordinates and sequences and names of substructures all the way to the right end utr_exon.";
}

sub process_breakpoint_list {
    my ($self, $breakpoints_list) = @_;
    my $out_fh = Genome::Sys->open_file_for_writing($self->fusion_output_file);
    die $self->error_message("fusion out file handle is missing") unless $out_fh;
    
    my @contents = map{'N/A'}column_names();
    my %output;

    for my $chr (keys %$breakpoints_list) {
        for my $item (@{$breakpoints_list->{$chr}}) {
            next if exists $item->{breakpoint_link};  
            my $aTranscriptRef = $item->{transcripts_crossing_breakpoint_a};
            my $bTranscriptRef = $item->{transcripts_crossing_breakpoint_b};

            my @item_types = qw(chrA bpA chrB bpB event);
            my $item_key   = $self->get_key_from_item($item);
            $output{$item_key} = \@contents;

	        # If there are no transcripts crossing either breakpoint, then you are done
	        next unless $aTranscriptRef and @$aTranscriptRef > 0 and $bTranscriptRef and @$bTranscriptRef > 0;
	    
	        # See if there are different or same genes crossing the breakpoints
            my $fusionProteinRef = $self->getAllInFrameFusions($aTranscriptRef, $item->{bpA}, $bTranscriptRef, $item->{bpB}, $item->{orient}, $out_fh);

	        if ( defined $fusionProteinRef && scalar(keys%{$fusionProteinRef}) >= 1 ) {
                my $aGeneNameRef = geneNames($aTranscriptRef);
	            my $bGeneNameRef = geneNames($bTranscriptRef);
                my $aGenes = join ' ', sort keys %$aGeneNameRef;
                my $bGenes = join ' ', sort keys %$bGeneNameRef;
                my $fusion = '['.$aGenes.']|['.$bGenes.']';

                my $line = join "\t", map{$item->{$_}}@item_types; 
                $out_fh->print("\n$line\t$fusion\n");
                $output{$item_key} = [$fusion];

		        for my $seq (keys %$fusionProteinRef) { 
		            for my $mRNA ( keys %{$fusionProteinRef->{$seq}} ) {
			            for my $mRNA_withCoordinates ( keys %{$fusionProteinRef->{$seq}->{$mRNA}} ) {
                            #$out_fh->print("\t$seq\n\t$mRNA\n\t$mRNA_withCoordinates\n"); 
                            $out_fh->print("\t$seq\n\t$mRNA_withCoordinates\n");
                        }
                    }
                }
            }
        }
    }
    $out_fh->close;
    return \%output;
}


sub getAllInFrameFusions {
    # Do all pairwise combinations of transcripts
    # Require both breakpoints to be in an intron (for now)
    # Check to see if relative orientation of sequences and genes make a fusion on one strand
    # Only include the fusions that are still in frame after putting the exons together

    my ($self, $aTranscriptRef, $bpA, $bTranscriptRef, $bpB, $junctionOrientation, $out_fh) = @_;
    my (%allFusions, %uniq); 
	  
    for my $aTranscript ( @$aTranscriptRef ) {
	    # Confirm breakpoint is in an intron
	    my $substructureA = $self->substructureWithBreakpoint($aTranscript, $bpA);
        next unless $substructureA and $substructureA =~ /intron/;

        my $a_gene_name       = $aTranscript->gene_name;
        my $a_transcript_name = $aTranscript->transcript_name;

	    for my $bTranscript ( @$bTranscriptRef ) {
	        # Confirm breakpoint is in an intron
	        my $substructureB = $self->substructureWithBreakpoint($bTranscript, $bpB);
            next unless $substructureB and $substructureB =~ /intron/;

            my $b_gene_name       = $bTranscript->gene_name;
            my $b_transcript_name = $bTranscript->transcript_name;
	        
            if ($a_gene_name eq $b_gene_name) {
                next unless $b_transcript_name eq $a_transcript_name; #Make sure same gene fusion wthin the same transcipt.
            }

	        # Determine which transcript is 5' and which is 3' of fused message based on
	        # relative orientation of sequences and strand that each gene is on
	        my $firstSecond = firstSecond($junctionOrientation, $aTranscript, $bTranscript);
	        # If transcripts are not going in same direction, there is no in-frame fusion possible and '$firstSecond' is returned undef
	        next unless defined $firstSecond;

	        # If both substructures are exons, get assembly contig and write to a file
	        # Need to check to see if event is CTX or other event and pass correct file to be read
	        # (i.e. need both assembly files as input so this won't work if a lot of different patient samples
	        #  are in file.  Will need to separate out patients or do file of files for assembly *fasta files.
	        #  Will also need to code an identifier for each patient.  Gets too complicated.
            my ($mRNA_a, $mRNA_b, $rnaWithCoordinates_a, $rnaWithCoordinates_b, $exon_ordinal_a, $exon_ordinal_b);

	        if ( $firstSecond eq "AB" ) {
		        ($mRNA_a, $rnaWithCoordinates_a, $exon_ordinal_a) = processedMessage($aTranscript, $bpA, "start");
		        ($mRNA_b, $rnaWithCoordinates_b, $exon_ordinal_b) = processedMessage($bTranscript, $bpB, "end");
		        $rnaWithCoordinates_a = "<".$a_gene_name."_". $aTranscript->transcript_name."> ".$rnaWithCoordinates_a;
		        $rnaWithCoordinates_b = "<".$b_gene_name."_". $bTranscript->transcript_name."> ".$rnaWithCoordinates_b;
            } 
            elsif ( $firstSecond eq "BA" ) {
		        ($mRNA_a, $rnaWithCoordinates_a, $exon_ordinal_a) = processedMessage($bTranscript, $bpB, "start");
		        ($mRNA_b, $rnaWithCoordinates_b, $exon_ordinal_b) = processedMessage($aTranscript, $bpA, "end");
		        $rnaWithCoordinates_a = "<".$b_gene_name."_". $bTranscript->transcript_name."> ".$rnaWithCoordinates_a;
		        $rnaWithCoordinates_b = "<".$a_gene_name."_". $aTranscript->transcript_name."> ".$rnaWithCoordinates_b;
            }
	        else {
		        confess "Unexpected \$firstSecond: '$firstSecond'";
            }

            next unless $mRNA_a and $mRNA_b;  #Must have both to proceed

            if ($a_gene_name eq $b_gene_name) {
                unless ($exon_ordinal_a and $exon_ordinal_b and $exon_ordinal_a =~ /^\d+$/ and $exon_ordinal_b =~ /^\d+$/) {
                    next;
                }
                next unless abs($exon_ordinal_a - $exon_ordinal_b) > 1; #so this works on same gene fusion (after deletion)
            }

	        my $fusionMessage = $mRNA_a.$mRNA_b;
	        next if $fusionMessage eq "";
	        my $seqObj  = new_sequence($fusionMessage);
	        my $protObj = $seqObj->translate;
	        my $protSeq = $protObj->seq;

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
		
		        # need to get translation of each end, but a codon could be split... first half could have different frame than second half
		        my $a_halfSeq = new_sequence($mRNA_a); 
                my $a_halfProtObj = $a_halfSeq->translate; 
                my $a_halfProtSeq = $a_halfProtObj->seq;
                my $b_halfProtSeq = inFrameSequence($mRNA_b);

                my $id = $a_halfProtSeq.'-'.$b_halfProtSeq;
                unless ($uniq{$id}) {
                    $out_fh->print("\n\n\$mRNA_a \t $a_halfProtSeq \n");
                    $out_fh->print("\$mRNA_b \t $b_halfProtSeq \n\n");
                    $uniq{$id} = 1;
                }

                if ($protSeq =~ /\Q$a_halfProtSeq\E/) {
                    my $new_seq = $a_halfProtSeq . ' | ';
                    $protSeq =~ s/\Q$a_halfProtSeq\E/$new_seq/;
                }
                $allFusions{$protSeq}{$message}{$messageWithCoordinates} = 1;  # put both protein sequences here, separated by '|' or ....
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
    my $gene           = $transcript->gene_name; 
    my $transcriptName = $transcript->transcript_name;
    my @ss             = $transcript->ordered_sub_structures;

    my $mRNA = "";
    my $rnaWithCoordinates = "";
    my $reading = 0;
    my $includeUtrExons = 1;  # for 3'
    if ( $startOrEnd eq "start" ) { 
	    $includeUtrExons = 0;  # for 5'
	    $reading = 1; 
    }

    my ($left, $right, $first_exon_ordinal, $last_exon_ordinal);

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
        #next if $sub_structure->structure_type eq "intron" or $sub_structure->structure_type eq 'flank';
        next unless $sub_structure->structure_type =~ /^(cds|utr)_exon$/;  #set this for now fir in-frame fusion
	
    	# The message to be translated has cds_exon and if part of the 5' has utr_exon that are beginning of that transcript
	    if ( $sub_structure->structure_type eq "cds_exon" ) {	     
	        $mRNA .= $sub_structure->nucleotide_seq; 
	        # Only include utr_exons until the first cds_exon
	        $includeUtrExons = 0; 
            $first_exon_ordinal = $sub_structure->ordinal unless $first_exon_ordinal;
            $last_exon_ordinal  = $sub_structure->ordinal;
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

    my $exon_ordinal = $startOrEnd eq "start" ? $last_exon_ordinal : $first_exon_ordinal;
    return ($mRNA, $rnaWithCoordinates, $exon_ordinal);
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
    
    my ($self, $transcript, $position) = @_;
    my @subStructures = $transcript->ordered_sub_structures;

    for my $structure (@subStructures) {
	    if ( $position >= $structure->{structure_start} and $position <= $structure->{structure_stop} ) {
	        return $structure->{structure_type}.$structure->{cds_exons_before};
        }
    }
    $self->warning_message("$position is not found in gene: ".$transcript->gene_name.'  transcript: '.$transcript->transcript_name);
    return;
}


sub column_names {
    return qw(fusion_genes);
}


