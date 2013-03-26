package Genome::Model::Tools::Annotate::Sv::Transcripts;

use strict;
use warnings;
use Carp;
use Genome;

class Genome::Model::Tools::Annotate::Sv::Transcripts {
    is => 'Genome::Model::Tools::Annotate::Sv::Base',
    has_input => [
        print_flanking_genes => {
            is => 'Boolean',
            default => 0,
            doc => 'Print columns that describe what genes fall within the flanking regions of the sv breakpoints',
        },
        cancer_gene_list => {
            is  => 'Text',
            doc => 'A list of cancer genes so they can be identified in output',
            is_optional => 1,
        },
    ],
    has_optional_transient => [
        _cancer_gene_names => {
            type => 'HASH',
        },
    ],
};

sub process_breakpoint_list{
    my ($self, $breakpoints_list) = @_;
    my $gene_list = $self->cancer_gene_list;

    if ($gene_list and -s $gene_list) {
        my $fh = Genome::Sys->open_file_for_reading($gene_list) or die "Failed to open $gene_list";
        my %genes;
        while (my $name = $fh->getline) {
            next if $name =~ /^#/;
            chomp $name;
            my @items = split /\s+/, $name;
            $genes{$items[0]} = 1;
        }
        $self->_cancer_gene_names(\%genes);
        $fh->close;
    }

    my %output;
    for my $chr (keys %$breakpoints_list) {
        for my $item (@{$breakpoints_list->{$chr}}) {
            next if $item->{breakpoint_link};
            my ($key, $value) = $self->process_item($item);
            $output{$key} = $value;
        }
    }
    return \%output;
}

sub process_item {
    my ($self, $item) = @_;

    my $geneA = my $transcriptA = my $orientationA = my $subStructureA = 'N/A';
    my $geneB = my $transcriptB = my $orientationB = my $subStructureB = 'N/A';
    my $deletedGenes = "N/A";
    my ($inCommonRef, $deletedGeneHashRef);

    # Get list of any deleted genes
    if ($item->{event} eq "DEL") { 
        my $transcriptRef = $item->{transcripts_between_breakpoints};
        $deletedGeneHashRef = geneNames($transcriptRef);
        if ( defined $deletedGeneHashRef and scalar(keys %$deletedGeneHashRef) > 1 ) {
            $deletedGenes = join ',', keys %$deletedGeneHashRef;
        }
    }
    # All transcripts crossing each breakpoint
    my $aTranscriptRef = $item->{transcripts_crossing_breakpoint_a};
    my $bTranscriptRef = $item->{transcripts_crossing_breakpoint_b};

    # If there are transcripts crossing both breakpoints see if there are transcripts in common
    if (defined $aTranscriptRef and scalar(@$aTranscriptRef) >= 1 and  defined $bTranscriptRef and scalar(@$bTranscriptRef) >= 1) {
        $inCommonRef = allTranscriptsInCommon($aTranscriptRef, $bTranscriptRef);

        # If there are transcripts in common, pick one to annotate
        # Preferentially choose one that affects protein product
        if ( defined $inCommonRef and scalar(@$inCommonRef) >= 1 ) {
            my $transcript = $self->bestTranscriptInCommon($inCommonRef, $item->{bpA}, $item->{bpB});
            (defined $transcript) || die "did not get transcript in common";
            ($geneA, $transcriptA, $orientationA, $subStructureA) = $self->transcriptFeatures($transcript, $item->{bpA});
            ($geneB, $transcriptB, $orientationB, $subStructureB) = $self->transcriptFeatures($transcript, $item->{bpB});	    
        }
    }

    # There are no transcripts in common.  Pick a different transcript for each breakpoint 
    if (defined $aTranscriptRef and scalar(@$aTranscriptRef) >= 1 and (not defined $inCommonRef or not scalar(@$inCommonRef))) { 
        my $transcript = $self->chooseBestTranscript($aTranscriptRef, $item->{bpA});
        ($geneA, $transcriptA, $orientationA, $subStructureA) = $self->transcriptFeatures($transcript, $item->{bpA});
    }
    if (defined $bTranscriptRef and scalar(@$bTranscriptRef) >= 1 and (not defined $inCommonRef or not scalar(@$inCommonRef))) { 
        my $transcript = $self->chooseBestTranscript($bTranscriptRef, $item->{bpB});
        ($geneB, $transcriptB, $orientationB, $subStructureB) = $self->transcriptFeatures($transcript, $item->{bpB});
    }

    my $key = $self->get_key_from_item($item);
    my $value = [$geneA, $transcriptA, $orientationA, $subStructureA, $geneB, $transcriptB, $orientationB, $subStructureB, $deletedGenes];

    if ($self->print_flanking_genes) {
        my $flankingARef  = $item->{transcripts_flanking_breakpoint_a};
        my $flankingBRef  = $item->{transcripts_flanking_breakpoint_b};

        my $flankingListA = _get_flanking_list($flankingARef);
        my $flankingListB = _get_flanking_list($flankingBRef);
       
        push @$value, $flankingListA;
        push @$value, $flankingListB;
    }

    return ($key, $value);
}

sub _get_flanking_list {
    my $list = shift;
    my $flankingList;

    if ($list and scalar @$list) {
        my %geneDedup;

        for my $t (@$list) {
            next unless $t;
            my $gene_name = $t->gene_name;
            $geneDedup{$gene_name} = 1 if $gene_name;  #Do we need to investigate the ones without name ?
        }
        $flankingList = join ",", keys %geneDedup if %geneDedup;
    }
    $flankingList = 'N/A' unless $flankingList;
    return $flankingList
}

sub allTranscriptsInCommon {
    # Return all transcripts that are in common 
    my ($aTranscriptRef, $bTranscriptRef) = @_;
    my (%aNames, @inCommon, $transcript);

    # Get a hash of names in first list
    for $transcript (@$aTranscriptRef) {
        my $transcriptName = $transcript->transcript_name;
        if (defined $transcriptName and $transcriptName ne "") {$aNames{$transcriptName} = 1;}
    }

    # See which names in second list are also in first    
    for $transcript (@$bTranscriptRef) {
        my $transcriptName = $transcript->transcript_name;
        if (defined $transcriptName and $transcriptName ne "" and  defined $aNames{$transcriptName}) {
            push @inCommon, $transcript;
        }
    }

    return \@inCommon;
}

sub geneNames {
    # Return ref to hash of gene names for the list of transcripts
    my $transcriptRef = shift;
    my %geneNames;
    for my $transcript (@$transcriptRef) {
        my $name = $transcript->gene_name;
        if ($name) { 
            $geneNames{$name} = 1; 
        }
    }
    return \%geneNames;
}

sub breakpointsInSameIntron {
    # Return 1 if the two breakpoints are in the same intron of the given transcript
    my ($self, $transcript, $bpA, $bpB) = @_;
    my $substructureA = $self->substructureWithBreakpoint($transcript, $bpA);
    my $substructureB = $self->substructureWithBreakpoint($transcript, $bpB);

    return unless $substructureA and $substructureB;
    return 1 if $substructureA eq $substructureB and $substructureA =~ /(flank)|(intron)/i;
}

sub substructureWithBreakpoint {
    # Return substructure with breakpoint.  The name of substructure is the type
    # and the number of coding exons before the substructure (does not utr_exon as coding)
    # type.cds_exons_before
    my ($self, $transcript, $position) = @_;
    my @subStructures = $transcript->ordered_sub_structures;

    for my $structure (@subStructures) {
        if ($position >= $$structure{structure_start} and $position <= $$structure{structure_stop}) {
            return $$structure{structure_type}.$$structure{cds_exons_before};
        }
    }
    $self->warning_message("$position is not found in gene: ".$transcript->gene_name.'  transcript: '.$transcript->transcript_name);
    return; #$position is not in any substructure of this transcript
}

sub transcriptFeatures {
    my ($self, $transcript, $position) = @_;

    my $geneName = $transcript->gene_name;
    if ($geneName) {
        my $cancer_gene_names = $self->_cancer_gene_names;
        if ($cancer_gene_names) {
            $geneName .= '(cancer)' if $cancer_gene_names->{$geneName};
        }
    }
    $geneName = 'N/A' unless $geneName;

    my $transcriptName = $transcript->transcript_name || 'N/A';
    my $orientation    = $transcript->strand || 'N/A';
    my $subStructure   = $self->substructureWithBreakpoint($transcript, $position) || 'N/A';
    
    return ($geneName, $transcriptName, $orientation, $subStructure);
}

sub bestTranscriptInCommon {
    # Choose transcript based on following priority:
    # 1. Gene product is affected (i.e. both breakpoints NOT in same intron) with gene name (there are a few Ensembl transcripts without one)
    # 2. Gene product is affected without a gene name 
    # 3. Transcript with gene name 
    # 4. Anything
    my ( $self, $inCommonRef, $bpA, $bpB ) = @_;
    my ( @notInSameIntron, @sameIntronWithGeneName );

    for my $transcript ( @$inCommonRef ) {
        my $geneName = $transcript->gene_name;
        if ( !$self->breakpointsInSameIntron($transcript, $bpA, $bpB) ) {
            if ( defined $geneName and $geneName ne "" ) { return $transcript; } # This is # 1. 
            push @notInSameIntron, $transcript; # This is # 2.
        } 
        elsif ( defined $geneName and $geneName ne "" ) {
            push @sameIntronWithGeneName, $transcript; # This is # 3.
        }
    }

    # If we are here, there is no # 1.
    if ( @notInSameIntron and scalar(@notInSameIntron) >= 1 ) { return $notInSameIntron[0]; }  # This is # 2.
    if ( @sameIntronWithGeneName and scalar(@sameIntronWithGeneName) >= 1 ) { return $sameIntronWithGeneName[0]; }  # This is # 3.
    return $inCommonRef->[0];
}

sub chooseBestTranscript {
    # Choose transcript based on following priority:
    # 1. Breakpoint in exon and has defined gene
    # 2. Want breakpoint in an exon
    # 3. Want transcript with gene name (there are a few Ensembl transcripts without one)
    # 4. Will take any transcript

    my ($self, $transcriptRef, $position) = @_;
    my ( @hasExon, @hasGeneName );

    for my $transcript ( @$transcriptRef ) { 
        my $substructure = $self->substructureWithBreakpoint($transcript, $position);
        next unless $substructure;
        my $geneName = $transcript->gene_name;
        if ( $substructure =~ /exon/i and defined $geneName and $geneName ne "" ) {
            # 1. This one has breakpoint in exon and has defined gene
            # This is all we need
            return $transcript;
        } 
        elsif ( $substructure =~ /exon/i ) {
            push @hasExon, $transcript;
        } 
        elsif ( defined $geneName and $geneName ne "" ) {
            push @hasGeneName, $transcript;
        }
    }

    # If we are here, there is no transcript with both breakpoint in exon
    # and gene name.  Send transcript based on following priority
    if ( @hasExon and scalar(@hasExon) >= 1 ) {
        # 2. Want breakpoint in an exon
        return $hasExon[0];
    } 
    elsif ( @hasGeneName and scalar(@hasGeneName) >= 1 ) {
        # 3. Want transcript with gene name
        return $hasGeneName[0];
    } 
    else {
        # 4. Will take any transcript
        return $transcriptRef->[0];
    }
    die "Should not be here";
}

sub find_annotated_positions {
    my ($self, $positions, $annotation, $annot_length, $tag) = @_;

    for my $chr (keys %$positions) {
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
                        for my $position_item (@{$positions->{$chr}}) {
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


sub column_names {
    my $self  = shift;
    my @names = qw(
        geneA transcriptA orientationA subStructureA 
        geneB transcriptB orientationB subStructureB	
        deletedGenes
    );
    if ($self->print_flanking_genes) {
        push @names, "flankingGenesA";
        push @names, "flankingGenesB";
    }
    return @names;
}
