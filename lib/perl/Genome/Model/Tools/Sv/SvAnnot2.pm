package Genome::Model::Tools::Sv::SvAnnot2;

use strict;
use warnings;
use Carp;
use Bio::Perl;
use Genome;

class Genome::Model::Tools::Sv::SvAnnot2 {
    is => 'Genome::Model::Tools::Sv',
    has => [
        breakdancer_files => {
            type => 'String',
            doc => 'Input files in breakdancer format',
            is_many => 1,
            is_optional => 1,
        },
        squaredancer_files => {
            type => 'String',
            doc => 'Input files in squaredancer format',
            is_many => 1,
            is_optional => 1,
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
    my $self = shift;
    my $build = $self->annotation_build;
    unless ($build){
        die "build not defined";
    }


# Get BreakDancer line, get first and second breakpoint
# Send each breakpoint to a sub that returns transcript and substructure containing the breakpoint
# What to do about using only one transcript for a given event.  Need to process them at same time (?)
# Just pick one transcript for a given gene.  And pick the same transcript for a given event.
# Get list of transcripts that cross a given breakpoint (in a sub)
# Then send both lists to another sub that picks one gene and one transcript in common if possible
# If not and both have gene, then it is fusion (?)
# But if DEL, then need everything between the transcripts

# Need to get affected trancript(s) for each bp

# 8.14    8       35255524        35255524        8       54284949        54284949        "DEL"


# all transcripts crossing each breakpoint
# transcripts in common
# foreach transcriptInCommon ==> substructure crossing breakpoint
#    if both substructures are same intron ==> no effect; otherwise ==> effect
#    choose any affected transcript with a gene name
    #

# affected: Both breakpoints NOT in same intron of a given transcript
#           Genes between breakpoints
    my $outfile = $self->output_file;
    open(OUT, "> $outfile") || die "Could not open output file: $!";
    print OUT "chrA\tbpA\tchrB\tbpB\tevent\tsize\tscore\tsource\tgeneA\ttranscriptA\torientationA\tsubStructureA\tgeneB\ttranscriptB\torientationB\tsubStructureB\tdeletedGenes\n";

    my @infiles = $self->breakdancer_files;
    for my $filename (@infiles) {
        processFile("bd", $filename, $build);
    }
    @infiles = $self->squaredancer_files;
    for my $filename (@infiles) {
        processFile("sd", $filename, $build);
    }

    close OUT;
}

sub processFile {
    my ($source,$infile, $build) = @_;

    my ( @entireFile, $line, $bdRef, $chrA, $bpA, $chrB, $bpB, $event, $patient, $aTranscriptRef, $bTranscriptRef, $transcriptRef,
        $eventId, $geneA, $transcriptA, $orientationA, $subStructureA, $geneB, $transcriptB, $orientationB, $subStructureB,
        $deletedGeneHashRef, $deletedGenes, $inCommonRef, $transcript, 
    );

    open(IN, "< $infile" ) || die "Could not open input file $infile: $!";
    @entireFile = <IN>;

    foreach my $line ( @entireFile ) {
        chomp $line;
        if ( $line =~ /^\s*$/ || $line =~ /^#/ ) { next; }
        if ( $line =~ /-------------------/) { last; } #to accomodate some of our squaredancer files with LSF output in them

        #parse either sd or bd file
        my ($id,$chrA,$bpA,$chrB,$bpB,$event,$oriA,$oriB,$size,$samples,$score);
        if ($source eq "bd") {
            ($id,$chrA,$bpA,undef,$chrB,$bpB,undef,$event,$oriA,undef,$size,$samples,$score) = split /\t/,$line;
            next if ($samples =~ /normal/);
        }
        elsif ($source eq "sd") {
            ($chrA,$bpA,undef,$chrB,$bpB,undef,$event,$size,$score) = split /\t/,$line;
        }

        my $index = 0;
        foreach my $var ($chrA,$bpA,$chrB,$bpB,$event,$size,$score) {
            unless (defined $var) { die "DID not define necessary variables for call:\n$line\n"; }
            $index++;
        }

        #$bdRef = BreakDancerLine->new($line);
        #($chrA, $bpA, $chrB, $bpB) = $bdRef->chromosomesAndBreakpoints();
        #$event = $bdRef->eventType(); 
        #$eventId = $bdRef->Id();
        $geneA = $transcriptA = $orientationA = $subStructureA = $geneB = $transcriptB = $orientationB = $subStructureB = $deletedGenes = "N/A";
        $inCommonRef = $deletedGeneHashRef = undef;

        #grab Ken's annotation
        #my @cols = split /\t/, $line;
        #my $kenAnnotation = $cols[14];

        # Get list of any deleted genes
        if ( $event eq "DEL" ) { 
            $transcriptRef = transcriptsBetweenBreakpoints($chrA, $bpA, $bpB, $build);
            $deletedGeneHashRef = geneNames($transcriptRef);
            if ( defined $deletedGeneHashRef && scalar(keys %{$deletedGeneHashRef}) > 1 ) {
                $deletedGenes = "";
                foreach (keys %{$deletedGeneHashRef}) { $deletedGenes .= "$_ "; }
            }
        }

        # All transcripts crossing each breakpoint
        $aTranscriptRef = allTranscripts($chrA, $bpA, $build);
        $bTranscriptRef = allTranscripts($chrB, $bpB, $build);

        # If there are no transcripts crossing the breakpoint, then you are done
        if ( (!defined $aTranscriptRef || scalar(@{$aTranscriptRef}) == 0) &&
            (!defined $bTranscriptRef || scalar(@{$bTranscriptRef}) == 0) ) {
            #print OUT "chrA\tbpA\tchrB\tbpB\tevent\tsize\tscore\tsource\tgeneA\ttranscriptA\torientationA\tsubStructureA\tgeneB\ttranscriptB\torientationB\tsubStructureB\tdeletedGenes\n";
            print OUT "$chrA\t$bpA\t$chrB\t$bpB\t$event\t$size\t$score\t$source\t$geneA\t$transcriptA\t$orientationA\t$subStructureA\t$geneB\t$transcriptB\t$orientationB\t$subStructureB\t$deletedGenes\n"; 
            next;
        }

        # If there are transcripts crossing both breakpoints see if there are transcripts in common
        if ( defined $aTranscriptRef && scalar(@{$aTranscriptRef}) >= 1 &&  defined $bTranscriptRef && scalar(@{$bTranscriptRef}) >= 1 ) {
            $inCommonRef = allTranscriptsInCommon($aTranscriptRef, $bTranscriptRef);
        }

        # If there are transcripts in common, pick one to annotate
        # Preferentially choose one that affects protein product
        if ( defined $inCommonRef && scalar(@{$inCommonRef}) >= 1 ) {
            $transcript = bestTranscriptInCommon($inCommonRef, $bpA, $bpB );
            (defined $transcript) || die "did not get transcript in common for \n'$line'";
            ($geneA, $transcriptA, $orientationA, $subStructureA) = transcriptFeatures($transcript, $bpA);
            ($geneB, $transcriptB, $orientationB, $subStructureB) = transcriptFeatures($transcript, $bpB);	    
            print OUT "$chrA\t$bpA\t$chrB\t$bpB\t$event\t$size\t$score\t$source\t$geneA\t$transcriptA\t$orientationA\t$subStructureA\t$geneB\t$transcriptB\t$orientationB\t$subStructureB\t$deletedGenes\n"; 
            next;	    
        }

        # There are no transcripts in common.  Pick a different transcript for each breakpoint 
        if ( defined $aTranscriptRef && scalar(@{$aTranscriptRef}) >= 1 ) { 
            $transcript = chooseBestTranscript($aTranscriptRef, $bpA);
            ($geneA, $transcriptA, $orientationA, $subStructureA) = transcriptFeatures($transcript, $bpA);
        }
        if ( defined $bTranscriptRef && scalar(@{$bTranscriptRef}) >= 1 ) { 
            $transcript = chooseBestTranscript($bTranscriptRef, $bpB);
            ($geneB, $transcriptB, $orientationB, $subStructureB) = transcriptFeatures($transcript, $bpB);
        }	    
        print OUT "$chrA\t$bpA\t$chrB\t$bpB\t$event\t$size\t$score\t$source\t$geneA\t$transcriptA\t$orientationA\t$subStructureA\t$geneB\t$transcriptB\t$orientationB\t$subStructureB\t$deletedGenes\n"; 

    } # match foreach $line

    close IN;

} # end of sub

sub allTranscripts {
    # Return all transcripts spanning the given position
    my ( $chr, $position, $build ) = @_;

    my ( $transcriptIterator, $transcript, @transcripts, );
    $transcriptIterator = $build->transcript_iterator(chrom_name => $chr); 
    (defined $transcriptIterator) || die "transcriptIterator not defined";
    while ( $transcript = $transcriptIterator->next ) {
        if ( $position >= $transcript->transcript_start() && $position <= $transcript->transcript_stop() ) {
            push @transcripts, $transcript;
        }
    }

    return \@transcripts;
}

sub allTranscriptsInCommon {
    # Return all transcripts that are in common 
    my ( $aTranscriptRef, $bTranscriptRef ) = @_;
    my ( %aNames, @inCommon, $transcript, $transcriptName,  );

    # Get a hash of names in first list
    foreach $transcript ( @{$aTranscriptRef} ) {
        $transcriptName = $transcript->transcript_name();
        if ( defined $transcriptName && $transcriptName ne "" ) {  $aNames{$transcriptName} = 1; }
    }

    # See which names in second list are also in first    
    foreach $transcript ( @{$bTranscriptRef} ) {
        $transcriptName = $transcript->transcript_name();
        if ( defined $transcriptName && $transcriptName ne "" &&  defined $aNames{$transcriptName} ) { push @inCommon, $transcript; }
    }

    return \@inCommon;
}


sub geneNames {
    # Return ref to hash of gene names for the list of transcripts
    my $transcriptRef = $_[0];
    my ( %geneNames, $transcript, $name,  );
    foreach $transcript ( @{$transcriptRef} ) {
        $name = $transcript->gene_name();
        if ( defined $name ) { $geneNames{$name} = 1; }
    }
    return \%geneNames;
}


sub breakpointsInSameIntron {
    # Return 1 if the two breakpoints are in the same intron of the given transcript
    my ( $transcript, $bpA, $bpB ) = @_;
    my $substructureA = substructureWithBreakpoint($transcript, $bpA);
    my $substructureB = substructureWithBreakpoint($transcript, $bpB);
    return ( $substructureA eq $substructureB && $substructureA =~ /(flank)|(intron)/i );
}

sub substructureWithBreakpoint {
    # Return substructure with breakpoint.  The name of substructure is the type
    # and the number of coding exons before the substructure (does not utr_exon as coding)
    # type.cds_exons_before

    my ( $transcript, $position ) = @_;

    my ( @subStructures, $structure, $subStructureDescription  );
    @subStructures = $transcript->ordered_sub_structures();
    foreach $structure ( @subStructures ) {
        if ( $position >= $$structure{structure_start} && $position <= $$structure{structure_stop} ) {
            return $$structure{structure_type}.$$structure{cds_exons_before};
        }
    }

    confess "'$position' is not in transcript";

}

sub transcriptFeatures {
    # Input: transcript, breakpoint
    # Return: $geneName, $transcriptName, $orientation, $subStructure
    my ( $transcript, $position ) = @_;
    my ( $geneName, $transcriptName, $subStructure, $orientation, );
    $geneName = $transcript->gene_name();
    if ( !defined $geneName || $geneName eq "" ) { $geneName = ""; }
    $transcriptName = $transcript->transcript_name();
    $subStructure = substructureWithBreakpoint($transcript, $position);
    $orientation = $transcript->strand();
    return ($geneName, $transcriptName, $orientation, $subStructure);
}

sub bestTranscriptInCommon {
    # Choose transcript based on following priority:
    # 1. Gene product is affected (i.e. both breakpoints NOT in same intron) with gene name (there are a few Ensembl transcripts without one)
    # 2. Gene product is affected without a gene name 
    # 3. Transcript with gene name 
    # 4. Anything

    my ( $inCommonRef, $bpA, $bpB ) = @_;
    my ( $transcript, @notInSameIntron, @sameIntronWithGeneName, $geneName, );
    foreach $transcript ( @{$inCommonRef} ) {
        $geneName = $transcript->gene_name();
        if ( !breakpointsInSameIntron($transcript, $bpA, $bpB) ) {
            if ( defined $geneName && $geneName ne "" ) { return $transcript; } # This is # 1. 
            push @notInSameIntron, $transcript; # This is # 2.
        } elsif ( defined $geneName && $geneName ne "" ) {
            push @sameIntronWithGeneName, $transcript; # This is # 3.
        }
    }

    # If we are here, there is no # 1.
    if ( @notInSameIntron && scalar(@notInSameIntron) >= 1 ) { return $notInSameIntron[0]; }  # This is # 2.
    if ( @sameIntronWithGeneName && scalar(@sameIntronWithGeneName) >= 1 ) { return $sameIntronWithGeneName[0]; }  # This is # 3.
    return $$inCommonRef[0];
}

sub chooseBestTranscript {
    # Choose transcript based on following priority:
    # 1. Breakpoint in exon and has defined gene
    # 2. Want breakpoint in an exon
    # 3. Want transcript with gene name (there are a few Ensembl transcripts without one)
    # 4. Will take any transcript

    my ($transcriptRef, $position) = @_;

    my ( $transcript, @hasExon, @hasGeneName, $geneName, $substructure,  );
    foreach $transcript ( @{$transcriptRef} ) { 
        $substructure = substructureWithBreakpoint($transcript, $position);
        $geneName = $transcript->gene_name();
        if ( $substructure =~ /exon/i && defined $geneName && $geneName ne "" ) {
            # 1. This one has breakpoint in exon and has defined gene
            #    This is all we need
            return $transcript;
        } elsif ( $substructure =~ /exon/i ) {
            push @hasExon, $transcript;
        } elsif ( defined $geneName && $geneName ne "" ) {
            push @hasGeneName, $transcript;
        }
    }

    # If we are here, there is no transcript with both breakpoint in exon
    # and gene name.  Send transcript based on following priority
    if ( @hasExon && scalar(@hasExon) >= 1 ) {
        # 2. Want breakpoint in an exon
        return $hasExon[0];
    } elsif ( @hasGeneName && scalar(@hasGeneName) >= 1 ) {
        # 3. Want transcript with gene name
        return $hasGeneName[0];
    } else {
        # 4. Will take any transcript
        return ${$transcriptRef}[0];
    }

    die "Should not be here";
}


sub transcriptsBetweenBreakpoints {
    # Used to get genes that are flanked by deletion breakpoints
    # The entire transcript has to be within the two breakpoints
    # Annotation of the individual breakpoints will give the transcripts interrupted by breakpoints
    my ($chr, $start, $stop, $build) = @_;
    if ( $start > $stop ) { ($start, $stop) = ($stop, $start); }
    my ( $transcriptIterator, $transcript, @transcripts, );

    $transcriptIterator = $build->transcript_iterator(chrom_name => $chr); 
    (defined $transcriptIterator) || die "transcriptIterator not defined";

    while ( $transcript = $transcriptIterator->next ) {
        #if ( $transcript->transcript_start() <= $stop && $transcript->transcript_stop() >= $start ) {

        if ( $transcript->transcript_start() >= $start && $transcript->transcript_start() <= $stop &&
            $transcript->transcript_stop()  >= $start && $transcript->transcript_stop() <= $stop ) {
            push @transcripts, $transcript;
        }
    }
    return \@transcripts;
}


sub processedMessage {
    # Return message up to given position
    my ( $transcript, $position ) = @_;
}
