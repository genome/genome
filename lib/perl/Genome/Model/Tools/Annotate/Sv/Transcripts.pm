package Genome::Model::Tools::Annotate::Sv::Transcripts;

use strict;
use warnings;
use Carp;
use Bio::Perl;
use Genome;

class Genome::Model::Tools::Annotate::Sv::Transcripts {
    is => 'Genome::Model::Tools::Annotate::Sv::Base',
    has => [
        breakpoint_wiggle_room => {
            is => 'Number',
            doc => 'Distance between breakpoint and annotated breakpoint within which they are considered the same, in bp',
            default => 200,
        },
#        dbsnp_annotation_file => {
#            is => 'Text',
#            doc => 'File containing UCSC dbsnp table',
#            default => "/gsc/scripts/share/BreakAnnot_file/human_build37/dbsnp132.indel.named.csv",
#        },
    ],
};

sub process_breakpoint_list{
    my ($self, $breakpoints_list) = @_;

    my %output;
    foreach my $chr (keys %{$breakpoints_list}) {
        foreach my $item (@{$breakpoints_list->{$chr}}) {
            my ($key, $value) = $self->process_item($item);
            $output{$key} = $value;
        }
    }
    return \%output;
    #my $dbsnp_annotation = $self->read_ucsc_annotation($self->dbsnp_annotation_file);

    #$self->find_annotated_positions($breakpoints_list, $dbsnp_annotation, $self->breakpoint_wiggle_room, "dbsnp_annotation");    
    
}

sub process_item {
    my $self = shift;
    my $item = shift;

    #only need to process each set of breakpoints once.
    if ($item->{breakpoint_link}) {
        return 1;
    }

    my ($geneA, $transcriptA, $orientationA, $subStructureA, $geneB, $transcriptB, $orientationB, $subStructureB, $deletedGenes, $aTranscriptRef, $bTranscriptRef, $inCommonRef, $deletedGeneHashRef);
    $geneA = $transcriptA = $orientationA = $subStructureA = $geneB = $transcriptB = $orientationB = $subStructureB = $deletedGenes = "N/A";
    $inCommonRef = $deletedGeneHashRef = undef;

    # Get list of any deleted genes
    if ( $item->{event} eq "DEL" ) { 
        my $transcriptRef = $item->{transcripts_between_breakpoints};
        $deletedGeneHashRef = geneNames($transcriptRef);
        if ( defined $deletedGeneHashRef and scalar(keys %{$deletedGeneHashRef}) > 1 ) {
            $deletedGenes = "";
            foreach (keys %{$deletedGeneHashRef}) { $deletedGenes .= "$_ "; }
        }
    }
    # All transcripts crossing each breakpoint
    $aTranscriptRef = $item->{transcripts_crossing_breakpoint_a};
    $bTranscriptRef = $item->{transcripts_crossing_breakpoint_b};

    # If there are transcripts crossing both breakpoints see if there are transcripts in common
    if ( defined $aTranscriptRef and scalar(@{$aTranscriptRef}) >= 1 and  defined $bTranscriptRef and scalar(@{$bTranscriptRef}) >= 1 ) {
        $inCommonRef = allTranscriptsInCommon($aTranscriptRef, $bTranscriptRef);

        # If there are transcripts in common, pick one to annotate
        # Preferentially choose one that affects protein product
        if ( defined $inCommonRef and scalar(@{$inCommonRef}) >= 1 ) {
            my $transcript = bestTranscriptInCommon($inCommonRef, $item->{bpA}, $item->{bpB} );
            (defined $transcript) || die "did not get transcript in common";
            ($geneA, $transcriptA, $orientationA, $subStructureA) = transcriptFeatures($transcript, $item->{bpA});
            ($geneB, $transcriptB, $orientationB, $subStructureB) = transcriptFeatures($transcript, $item->{bpB});	    
        }
    }

    # There are no transcripts in common.  Pick a different transcript for each breakpoint 
    if ( defined $aTranscriptRef and scalar(@{$aTranscriptRef}) >= 1 and (not defined $inCommonRef or not scalar(@$inCommonRef))) { 
        my $transcript = chooseBestTranscript($aTranscriptRef, $item->{bpA});
        ($geneA, $transcriptA, $orientationA, $subStructureA) = transcriptFeatures($transcript, $item->{bpA});
    }
    if ( defined $bTranscriptRef and scalar(@{$bTranscriptRef}) >= 1 and (not defined $inCommonRef or not scalar(@$inCommonRef))) { 
        my $transcript = chooseBestTranscript($bTranscriptRef, $item->{bpB});
        ($geneB, $transcriptB, $orientationB, $subStructureB) = transcriptFeatures($transcript, $item->{bpB});
    }

    my $dbsnp_ref = $item->{dbsnp_annotation};
    my $dbsnp_string = "-";
    if ($dbsnp_ref) {
        my @dbsnp = map {$_->{name}} @{$dbsnp_ref};
        $dbsnp_string = join(",", @dbsnp);
    }
    my $key = join("--", $item->{chrA}, $item->{bpA}, $item->{chrB}, $item->{bpB}, $item->{event});
    my $value = [$geneA, $transcriptA, $orientationA, $subStructureA, $geneB, $transcriptB, $orientationB, $subStructureB, $deletedGenes];
    return ($key, $value);

    return 1;
}

sub allTranscriptsInCommon {
    # Return all transcripts that are in common 
    my ( $aTranscriptRef, $bTranscriptRef ) = @_;
    my ( %aNames, @inCommon, $transcript, $transcriptName,  );

    # Get a hash of names in first list
    foreach $transcript ( @{$aTranscriptRef} ) {
        $transcriptName = $transcript->transcript_name();
        if ( defined $transcriptName and $transcriptName ne "" ) {  $aNames{$transcriptName} = 1; }
    }

    # See which names in second list are also in first    
    foreach $transcript ( @{$bTranscriptRef} ) {
        $transcriptName = $transcript->transcript_name();
        if ( defined $transcriptName and $transcriptName ne "" and  defined $aNames{$transcriptName} ) { push @inCommon, $transcript; }
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
    return ( $substructureA eq $substructureB and $substructureA =~ /(flank)|(intron)/i );
}

sub substructureWithBreakpoint {
    # Return substructure with breakpoint.  The name of substructure is the type
    # and the number of coding exons before the substructure (does not utr_exon as coding)
    # type.cds_exons_before

    my ( $transcript, $position ) = @_;

    my ( @subStructures, $structure, $subStructureDescription  );
    @subStructures = $transcript->ordered_sub_structures();
    foreach $structure ( @subStructures ) {
        if ( $position >= $$structure{structure_start} and $position <= $$structure{structure_stop} ) {
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
            if ( defined $geneName and $geneName ne "" ) { return $transcript; } # This is # 1. 
            push @notInSameIntron, $transcript; # This is # 2.
        } elsif ( defined $geneName and $geneName ne "" ) {
            push @sameIntronWithGeneName, $transcript; # This is # 3.
        }
    }

    # If we are here, there is no # 1.
    if ( @notInSameIntron and scalar(@notInSameIntron) >= 1 ) { return $notInSameIntron[0]; }  # This is # 2.
    if ( @sameIntronWithGeneName and scalar(@sameIntronWithGeneName) >= 1 ) { return $sameIntronWithGeneName[0]; }  # This is # 3.
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
        if ( $substructure =~ /exon/i and defined $geneName and $geneName ne "" ) {
            # 1. This one has breakpoint in exon and has defined gene
            #    This is all we need
            return $transcript;
        } elsif ( $substructure =~ /exon/i ) {
            push @hasExon, $transcript;
        } elsif ( defined $geneName and $geneName ne "" ) {
            push @hasGeneName, $transcript;
        }
    }

    # If we are here, there is no transcript with both breakpoint in exon
    # and gene name.  Send transcript based on following priority
    if ( @hasExon and scalar(@hasExon) >= 1 ) {
        # 2. Want breakpoint in an exon
        return $hasExon[0];
    } elsif ( @hasGeneName and scalar(@hasGeneName) >= 1 ) {
        # 3. Want transcript with gene name
        return $hasGeneName[0];
    } else {
        # 4. Will take any transcript
        return ${$transcriptRef}[0];
    }

    die "Should not be here";
}

sub find_annotated_positions {
    my $self = shift;
    my $positions = shift;
    my $annotation = shift;
    my $annot_length = shift;
    my $tag = shift;

    foreach my $chr (keys %$positions) {
        my @sorted_items = sort {$a->{bpB}<=>$b->{bpB}} (@{$positions->{$chr}});
        my @sorted_positions = map{$_->{bpB}} @sorted_items;
        my %annotated_output;

        my @chromEnds = sort {$a<=>$b} keys %{$annotation->{$chr}};
        for my $pos (@sorted_positions) {
            while (@chromEnds>0 && $pos>$chromEnds[0]+$annot_length) {
                shift @chromEnds;
            }
            next unless @chromEnds>0;
            for my $start (keys %{$$annotation{$chr}{$chromEnds[0]}}) {
                if ($pos>=$start-$annot_length) {
                    for my $var (@{$$annotation{$chr}{$chromEnds[0]}{$start}}){
                        foreach my $position_item (@{$positions->{$chr}}) {
                            if ($position_item->{bpB} eq $pos) {
                                push @{$position_item->{$tag}}, $var;
                            }
                        }
                    }
                }
            }
        }
    }
    return 1;
}

sub read_ucsc_annotation{
    my ($self, $file) = @_;
    my %annotation;
    open (ANNOTATION, "<$file") || die "Unable to open annotation: $file\n";
    while (<ANNOTATION>) {
        chomp;
        next if /^\#/;
        my $p;
        my @extra;
        ($p->{bin},$p->{chrom},$p->{chromStart},$p->{chromEnd},$p->{name},@extra) = split /\t+/;
        $p->{chrom} =~ s/chr//;
        $p->{extra} = \@extra;
        push @{$annotation{$p->{chrom}}{$p->{chromEnd}}{$p->{chromStart}}}, $p;
    }
    close ANNOTATION;
    return \%annotation;
}


sub column_names {
    return qw(
        geneA transcriptA orientationA subStructureA 
        geneB transcriptB orientationB subStructureB	
        deletedGenes
    );
}
