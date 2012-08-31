package Genome::Model::Tools::Sv::AssemblyPipeline::Hits;

use Genome;

class Genome::Model::Tools::Sv::AssemblyPipeline::Hits {
    is => 'Genome::Model::Tools::Sv::AssemblyPipeline'
};

use strict;
use warnings;
use Carp;


# This is used to parse the cross_match hits
my $number = "\\d+\\.?\\d*";
my $deleted = "\\(\\s*\\d+\\s*\\)";
my $AlignmentLine = "($number)\\s+($number)\\s+($number)\\s+($number)\\s+(\\S+)\\s+(\\d+)\\s+(\\d+)\\s+($deleted)\\s+(C\\s+)?(\\S+)\\s+($deleted)?\\s*(\\d+)\\s+(\\d+)\\s*($deleted)?";


########
#
#
#   need to deal with output from '-discrep_lists' flag
#   
#





sub new {
    my $proto = shift;
    my $class = ref($proto) || $proto;
    my $self = {};

    $self->{ALIGNMENT} = undef; 
    $self->{SCORE} = undef; 
    $self->{SUBS} = undef;
    $self->{DELS} = undef;
    $self->{INS} = undef;
    $self->{QUERY} = undef;
    $self->{QSTART} = undef;
    $self->{QEND} = undef;
    $self->{PAST_QUERY_END} = undef;
    $self->{COMPLEMENT} = 0;
    $self->{SUBJECT} = undef;
    $self->{BEFORE_SUBJECT_START} = undef;
    $self->{S_START} = undef;
    $self->{S_END} = undef;
    $self->{AFTER_SUBJECT_END} = undef;
    $self->{NOT_BEST_HIT} = 0;      #  there is a higher-scoring match whose domain partly 
                                    #  includes the domain of this match

    $self->{INSERTIONS} = ();  # Hash of insertions in alignment. ${$hash{$queryStart}}{$targetStart} = $length
    $self->{DELETIONS} = ();   # Hash of deletions in alignment. ${$hash{$queryStart}}{$targetStart} = $length

    bless ($self, $class);
    return $self;
}

sub addCrossMatchLine {
    my $self = shift;
    my $line = shift;

    if ( $line =~ /$AlignmentLine/ ) {
	$self->score($1);
	$self->percentSubs($2);
	$self->percentDels($3);
	$self->percentInserts($4);
	$self->queryName($5);
	$self->queryStart($6);
	$self->queryEnd($7);
	$self->pastQueryEnd($8);
	$self->isComplemented($9);
	$self->subjectName($10);
	$self->beforeSubjectStart($11);
	$self->subjectStart($12);
	$self->subjectEnd($13);
	$self->afterSubjectEnd($14);
	$self->alignmentLine($line);
    } else {
	confess "Unexpected format for alignment line: '$line'";
    }

    if ( $line =~ /\*/ ) { $self->notBestHit(1); }
}

sub alignmentLine {
    my $self = shift;
    if (@_) { $self->{ALIGNMENT} = shift; }
    return $self->{ALIGNMENT};
}

sub score {
    my $self = shift;
    if (@_) { $self->{SCORE} = shift; }
    return $self->{SCORE};
}

sub percentSubs {
    my $self = shift;
    if (@_) { $self->{SUBS} = shift; }
    return $self->{SUBS};
}

sub percentDels {
    my $self = shift;
    if (@_) { $self->{DELS} = shift; }
    return $self->{DELS};
}

sub percentInserts {
    my $self = shift;
    if (@_) { $self->{INS} = shift; }
    return $self->{INS};
}

sub queryName {
    my $self = shift;
    if (@_) { $self->{QUERY} = shift; }
    return $self->{QUERY};
}

sub queryStart {
    my $self = shift;
    if (@_) { $self->{QSTART} = shift; }
    return $self->{QSTART};
}

sub queryEnd {
    my $self = shift;
    if (@_) { $self->{QEND} = shift; }
    return $self->{QEND};
}

sub pastQueryEnd {
    my $self = shift;
    if (@_) { 
	$self->{PAST_QUERY_END} = shift; 
	if (defined $self->{PAST_QUERY_END} && $self->{PAST_QUERY_END} =~ /(\d+)/ ) { 
	    $self->{PAST_QUERY_END} = $1; 
	}
    }
    return $self->{PAST_QUERY_END};
}

sub querySize {
    my $self = shift;
    
    my $size = $self->pastQueryEnd() + $self->queryEnd();
    return $size;
}

sub subjectSize {
    my $self = shift;
    my $subjectSize;
    if ( $self->isComplemented() ) {
	$subjectSize = $self->subjectStart() + $self->beforeSubjectStart();
    } else {
	$subjectSize = $self->subjectEnd() + $self->afterSubjectEnd();
    }	  
    
    return $subjectSize;
}


sub isComplemented {
    my $self = shift;
    if (@_) { 
	$self->{COMPLEMENT} = shift; 
	if ( defined $self->{COMPLEMENT} && $self->{COMPLEMENT} =~ /C/ ) {  
	    $self->{COMPLEMENT} = 1; 
	} else { 
	    $self->{COMPLEMENT} = 0; 
	}  
    }
    return $self->{COMPLEMENT};
}


sub subjectName {
    my $self = shift;
    if (@_) { $self->{SUBJECT} = shift; }
    return $self->{SUBJECT};
}

sub beforeSubjectStart  {
    my $self = shift;
    if (@_) { 
	$self->{BEFORE_SUBJECT_START} = shift; 
	if ( defined $self->{BEFORE_SUBJECT_START} && $self->{BEFORE_SUBJECT_START} =~ /(\d+)/ ) { 
	    $self->{BEFORE_SUBJECT_START} = $1; 
	}
    }
    return $self->{BEFORE_SUBJECT_START};
}

sub subjectStart {
    my $self = shift;
    if (@_) { $self->{S_START} = shift; }
    return $self->{S_START};
}

sub subjectEnd {
    my $self = shift;
    if (@_) { $self->{S_END} = shift; }
    return $self->{S_END};
}

sub afterSubjectEnd {
    my $self = shift;
    if (@_) { 
	$self->{AFTER_SUBJECT_END} = shift; 
	if ( defined $self->{AFTER_SUBJECT_END} && $self->{AFTER_SUBJECT_END} =~ /(\d+)/ ) {  
	    $self->{AFTER_SUBJECT_END} = $1; 
	}
    }
    return $self->{AFTER_SUBJECT_END};
}

sub alignmentsAreSame {
    # See if two alignments are the same (except for query name)
    # Used to de-duplicate reads (reads with same sequence will have same alignment to same target)
    my ($self, $hitRef) = @_;

    # Get the alignment line for $self and $hitRef.  Everything should be the same except for the 
    # query name (which is in the 5th column; 4th using first as 0)

    my $selfAlignmentLine = $self->alignmentLine();
    my $hitRefAlignmentLine = $hitRef->alignmentLine();

    # Get rid of any white space at beginning or end
    $hitRefAlignmentLine =~ s/^\s*//; $hitRefAlignmentLine =~ s/\s*$//;
    $selfAlignmentLine =~ s/^\s*//; $selfAlignmentLine =~ s/\s*$//;

    my @selfColumns = split /\s+/, $selfAlignmentLine;
    my @otherColumns = split /\s+/, $hitRefAlignmentLine;
   
    for ( my $i = 0; $i <= $#otherColumns; $i++ ) {
	if ( $i != 4 && $selfColumns[$i] ne $otherColumns[$i] ) { return 0; }
    }

    return 1;
}
 

sub notBestHit {
    my $self = shift;
    if (@_) { $self->{NOT_BEST_HIT} = shift; }
    return $self->{NOT_BEST_HIT};
}

sub addInsertion {
    # This is an insertion within a single alignment
    # Get this information using '-discrep_lists' flag in cross_match
    my $self = shift;
    my ($queryStart, $targetStart, $length) = @_;
    ( !defined ${${$self->{INSERTIONS}}{$queryStart}}{$targetStart} ) ||
	carp "Duplicate insertion entry at \$queryStart = $queryStart in $self->{ALIGNMENT}";
    
    ${${$self->{INSERTIONS}}{$queryStart}}{$targetStart} = $length;
}

sub addDeletion {
    # This is a deletion within a single alignment
    # Get this information using '-discrep_lists' flag in cross_match
    my $self = shift;
    my ($queryStart, $targetStart, $length) = @_;
    ( !defined ${${$self->{DELETIONS}}{$queryStart}}{$targetStart} ) ||
	carp "Duplicate deletion entry at \$queryStart = $queryStart in $self->{ALIGNMENT}";

    ${${$self->{DELETIONS}}{$queryStart}}{$targetStart} = $length;
}

sub returnDeletions {
    # return ref to hash of hash with keys $queryStart, $targetStart, value = length
    #     ${$$hash{$queryStart}}{$targetStart} = $length
    my $self = shift;
    return \%{$self->{DELETIONS}};
}

sub returnInsertions {
    # return ref to hash of hash with keys $queryStart, $targetStart, value = length
    #     ${$$hash{$queryStart}}{$targetStart} = $length
    my $self = shift;
    return \%{$self->{INSERTIONS}};
}
 
sub returnLargestDeletion {
    # Return ($delQueryStart, $delTargetStart, $largestLength) from the 
    # largest deletion (if there is one)

    my $self = shift;
    my $indelRef = $self->returnDeletions();
    if ( !defined $indelRef || scalar(keys %{$indelRef}) == 0 ) { return undef; }

    my ($queryStart, $targetStart,  $delQueryStart, $delTargetStart, $largestLength, %largestDel, );
    $largestLength = 0;
    
    foreach $queryStart (keys %{$indelRef}) {
	foreach $targetStart ( keys %{$$indelRef{$queryStart}} ) {
	    if ( ${$$indelRef{$queryStart}}{$targetStart} > $largestLength ) {
		$largestLength = ${$$indelRef{$queryStart}}{$targetStart};
		$delQueryStart = $queryStart;
		$delTargetStart = $targetStart;
	    }
	}
    }

    (defined $delQueryStart && defined $delTargetStart) ||
	confess "Expected del start and stop to be defined", $self->alignmentLine();

    return ($delQueryStart, $delTargetStart, $largestLength);
}

sub returnDeletionTotal {
    # Return the sum of all deletion lengths in this hit
    my $self = shift;
    my ($queryStart, $targetStart, $sum, );
    my $indelRef = $self->returnDeletions();
    if ( !defined $indelRef || scalar(keys %{$indelRef}) == 0 ) { return 0; }
    $sum = 0;
    foreach $queryStart (keys %{$indelRef}) {
	foreach $targetStart ( keys %{$$indelRef{$queryStart}} ) {
	    $sum += ${$$indelRef{$queryStart}}{$targetStart};
	}
    }
    
    return $sum;
}

sub returnInsertionTotal {
    # Return the sum of all insertion lengths in this hit
    my $self = shift;
    my ($queryStart, $targetStart, $sum, );
    my $indelRef = $self->returnInsertions();
    if ( !defined $indelRef || scalar(keys %{$indelRef}) == 0 ) { return 0; }
    $sum = 0;
    foreach $queryStart (keys %{$indelRef}) {
	foreach $targetStart ( keys %{$$indelRef{$queryStart}} ) {
	    $sum += ${$$indelRef{$queryStart}}{$targetStart};
	}
    }
    
    return $sum;
}


sub returnLargestInsertion {
    # Return ($insQueryStart, $insTargetStart, $largestLength) from the 
    # largest insertion (if there is one)

    my $self = shift;
    my $indelRef = $self->returnInsertions();
    if ( !defined $indelRef || scalar(keys %{$indelRef}) == 0 ) { return undef; }

    my ($queryStart, $targetStart,  $insQueryStart, $insTargetStart, $largestLength, %largestIns, );
    $largestLength = 0;
    
    foreach $queryStart (keys %{$indelRef}) {
	foreach $targetStart ( keys %{$$indelRef{$queryStart}} ) {
	    if ( ${$$indelRef{$queryStart}}{$targetStart} > $largestLength ) {
		$largestLength = ${$$indelRef{$queryStart}}{$targetStart};
		$insQueryStart = $queryStart;
		$insTargetStart = $targetStart;
	    }
	}
    }

    (defined $insQueryStart && defined $insTargetStart) ||
	confess "Expected ins start and stop to be defined", $self->alignmentLine();
   
    return ($insQueryStart, $insTargetStart, $largestLength);
}



sub deletionReferenceBoundaries {
    # Return: start and stop of deletion on target (i.e. min($targetStart), max($targetStart + $length) )
    # This is not ideal since there might be small outlier deletions.  Maybe better to get
    # the single largest indel?? (but sometimes there can be a superious short alignment in the middle)
    #   return ($minTargetStart, $maxTargetStart);

    my $self = shift;
    my $indelRef = $self->returnDeletions();
    my ( $minTargetStart, $maxTargetStart, );
    
    carp "Consider using Hits::returnLargestDeletion";

    $minTargetStart = $maxTargetStart = undef;
    foreach my $queryStart (keys %{$indelRef}) {
	foreach my $targetStart ( keys %{$$indelRef{$queryStart}} ) {
	    if ( !defined $minTargetStart || $targetStart < $minTargetStart ) {
		$minTargetStart = $targetStart;
	    }
	    if ( !defined $maxTargetStart || 
		 $targetStart + ${$$indelRef{$queryStart}}{$targetStart} >  $maxTargetStart ) {
		$maxTargetStart = $targetStart + ${$$indelRef{$queryStart}}{$targetStart};
	    }
	} # match foreach my $targetStart
    } # match foreach my $queryStart

    return ($minTargetStart, $maxTargetStart);
}

sub deletionContigBoundaries {
    # Return: start and stop of deletion on contig (i.e. min($queryStart), max($queryStart) )
    #         The start and stop should be about the same
    # This is not ideal since there might be small outlier deletions.  Maybe better to get
    # the single largest indel?? (but sometimes there can be a superious short alignment in the middle)
    #   return ($minQueryStart, $maxQueryStart);

    my $self = shift;
    my $indelRef = $self->returnDeletions();
    my ( $minQueryStart, $maxQueryStart, );
    
    carp "Consider using Hits::returnLargestDeletion";

    $minQueryStart = $maxQueryStart = undef;
    foreach my $queryStart (keys %{$indelRef}) {
	if ( !defined $minQueryStart || $queryStart < $minQueryStart ) {
	    $minQueryStart = $queryStart;
	}
	if ( !defined $maxQueryStart || $queryStart > $maxQueryStart ) {
	    $maxQueryStart = $queryStart;
	}
    } # match foreach my $queryStart

    return ($minQueryStart, $maxQueryStart);
}
    

sub insertionReferenceBoundaries {
    # Return: start and stop of insertion on target (i.e. min($targetStart), max($targetStart) )
    #         The start and stop should be about the same
    # This is not ideal since there might be small outlier deletions.  Maybe better to get
    # the single largest indel?? (but sometimes there can be a superious short alignment in the middle)
    #   return ($minTargetStart, $maxTargetStart);

    my $self = shift;
    my $indelRef = $self->returnInsertions;
    my ( $minTargetStart, $maxTargetStart, );
    

    
    carp "Consider using Hits::returnLargestInsertion";

    $minTargetStart = $maxTargetStart = undef;
    foreach my $queryStart (keys %{$indelRef}) {
	foreach my $targetStart ( keys %{$$indelRef{$queryStart}} ) {
	    if ( !defined $minTargetStart || $targetStart < $minTargetStart ) {
		$minTargetStart = $targetStart;
	    }
	    if ( !defined $maxTargetStart || 
		 $targetStart >  $maxTargetStart ) {
		$maxTargetStart = $targetStart;
	    }
	} # match foreach my $targetStart
    } # match foreach my $queryStart

    return ($minTargetStart, $maxTargetStart);
}
   
 

sub insertionContigBoundaries {
    # Return: start and stop of insertion on contig (i.e. min($queryStart), max($queryStart + $length) )
    # This is not ideal since there might be small outlier deletions.  Maybe better to get
    # the single largest indel?? (but sometimes there can be a superious short alignment in the middle)
    #   return ($minQueryStart, $maxQueryStart);

    my $self = shift;
    my $indelRef = $self->returnInsertions;
    my ( $minQueryStart, $maxQueryStart, );
    
    carp "Consider using Hits::returnLargestInsertion";

    $minQueryStart = $maxQueryStart = undef;
    foreach my $queryStart (keys %{$indelRef}) {
	foreach my $targetStart ( keys %{$$indelRef{$queryStart}} ) {
	    if ( !defined $minQueryStart || $queryStart < $minQueryStart ) {
		$minQueryStart = $queryStart;
	    }
	    if ( !defined $maxQueryStart || 
		 $queryStart + ${$$indelRef{$queryStart}}{$targetStart} >  $maxQueryStart ) {
		$maxQueryStart = $queryStart + ${$$indelRef{$queryStart}}{$targetStart};
	    }
	} # match foreach my $queryStart
    } # match foreach my $queryStart

    return ($minQueryStart, $maxQueryStart);
}
   





1;
