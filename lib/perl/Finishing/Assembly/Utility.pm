package Finishing::Assembly::Utility;

use Exporter;
use Compress::Zlib;
use List::Util qw(min max);
@ISA    = qw(Exporter);
@EXPORT = qw(
calculate_intersection
intervals_intersect
check_contiguous

);

our $VERSION = 0.01;

use strict;

#use warnings;
use Carp;
use Carp::Assert;
use Cwd;
use IO::String;

use constant QUALITY_MIN               => 0;
use constant QUALITY_LOW               => 0;
use constant QUALITY_WITHOUT_PHD_FILES => 15;
use constant DEF_QUAL_VALUE            => 20;
use constant QUALITY_HIGH              => 50;
use constant QUALITY_LOW_EDITED        => 98;
use constant QUALITY_HIGH_EDITED       => 99;
use constant DEFAULT_QUALITY_EDITED    => QUALITY_LOW_EDITED;
use constant QUALITY_MAX               => 100;

use constant TOP_STRAND_TERMINATOR    => 0;
use constant TOP_STRAND_PRIMER        => 1;
use constant BOTTOM_STRAND_TERMINATOR => 2;
use constant BOTTOM_STRAND_PRIMER     => 3;

#functions for padding and unpadding consensus and phd quality arrays

sub nNormalQualityFrom9899Quality
{
	my ($nQuality) = @_;
	if ( $nQuality == QUALITY_LOW_EDITED )
	{
		$nQuality = QUALITY_MIN;
	}
	if ( $nQuality == QUALITY_HIGH_EDITED )
	{
		$nQuality = 60;
	}

	# Phil wanted it this way Dec 1997
	# for data quality exercise to increase the
	# quality of pads between base 99's
	return $nQuality;
}

sub _bIsDyePrimerNotDyeTerminator
{
	my ($rRead) = @_;

	if ( defined( $rRead->{chem} ) )
	{
		if ( $rRead->{chem} =~ /prim/ )
		{
			return 1;
		}
		elsif ( $rRead->{chem} =~ /term/ )
		{
			return 0;
		}
	}

	# if reached here, there are 2 possibilities:
	# Either the CHEM field was not set
	# or, if it was, it was set to something other than
	# prim or term (such as unknown or something unrecognizeable )

	# Thus use the filename extension to detect it.
	my $sFileName = $rRead->{chromat_file};

	my ( $Name, $Ext ) = $sFileName =~ /^(.+\.)(.+)$/;

	if ( !defined($Ext) )
	{
		return 1;
	}

	if ( ( substr( $Ext, 0, 1 ) =~ /s|r|q/ ) )
	{
		return 1;
	}
	else
	{
		return 0;
	}
}

sub _nGetStrandAndChemistry
{
	my ($rRead) = @_;
	my $bPrimerNotTerminator = _bIsDyePrimerNotDyeTerminator($rRead);

	my $bBottomStrand = $rRead->{complemented};

	if ($bPrimerNotTerminator)
	{
		if ($bBottomStrand)
		{
			return (BOTTOM_STRAND_PRIMER);
		}
		else
		{
			return (TOP_STRAND_PRIMER);
		}
	}
	else
	{
		if ($bBottomStrand)
		{
			return (BOTTOM_STRAND_TERMINATOR);
		}
		else
		{
			return (TOP_STRAND_TERMINATOR);
		}
	}
}

sub _nFindIndexOfMatchOrPredecessor
{
	my ( $rArray, $nIndex ) = @_;

	if ( $rArray == 0 )
	{
		return undef;
	}

	# region A
	# ----------- nMinIndex
	# region B
	# ----------- nTestIndex
	# region C
	# ----------- nTooLargeIndex
	# region D

	# an index becomes nTooLargeIndex if it is greater than soMatch
	# an index becomes nMinIndex if it is less than soMatch

	return undef if $nIndex < $$rArray[0];

	my $nMinIndex    = 0;
	my $nTooBigIndex = @$rArray - 1;

	if (   $$rArray[$nTooBigIndex] < $nIndex
		|| $$rArray[$nTooBigIndex] == $nIndex )
	{
		return $nTooBigIndex;
	}

	# if reached here, soMatch < operator[]( nTooBigIndex )
	# and operator[](0) <= soMatch
	# Thus nTooBigIndex != 0  and the correct index is somewhere
	# less than nTooBigIndex and correct index >= 0 so
	# correct index >= 0.
	# Thus bisect this range over and over, making it smaller
	# and smaller until it is 0.

	while (1)
	{
		if ( $nTooBigIndex - $nMinIndex <= 1 )
		{
			return ($nMinIndex);
		}
		else
		{
			my $nTestIndex = ( $nTooBigIndex + $nMinIndex ) / 2;
			if ( $nIndex < $$rArray[$nTestIndex] )
			{
				$nTooBigIndex = $nTestIndex;
			}
			else
			{
				$nMinIndex = $nTestIndex;
			}
		}
	}
}

sub _bIsWholeCloneTemplateNotSubcloneTemplate
{
	my ($rRead) = @_;
	if ( defined( $rRead->{template} ) )    #check phd template first
	{

		if ( $rRead->{template} =~ /bac|cos/ )
		{
			return 1;
		}
		else
		{
			return 0;
		}
	}

	# if reached here, then the information was not in the phd file

	if ( $rRead->{name} =~ /_BAC|_bac/ )
	{
		return 1;
	}
	else
	{
		return 0;
	}
}

# so far, this is used to tell whether 2 reads come from the
# same subclone template, and thus may have the same mutation
sub _bAreTheseReadsFromTheSameTemplate
{

	my ( $rRead1, $rRead2 ) = @_;

	if (   _bIsWholeCloneTemplateNotSubcloneTemplate($rRead1)
		|| _bIsWholeCloneTemplateNotSubcloneTemplate($rRead2) )
	{
		return 0;
	}

	# if reached here, neither are whole clone templates

	if ( $rRead1->{template} eq $rRead2->{template} )
	{
		return 1;
	}
	else
	{
		return 0;
	}
}

sub _cons_position_is_in_high_quality_segment_of_read
{
	my ( $read, $cons_pos ) = @_;
	if ( $read->{qual_clip_start} == -1 && $read->{qual_clip_stop} == -1 )
	{
		return 0;
	}
	else
	{
		if ( $read->{qual_clip_start} <= get_child_position_from_parent_position($read, $cons_pos )
			&& get_child_position_from_parent_position($read, $cons_pos )<= $read->{qual_clip_stop} )
		{
			return 1;
		}
		else
		{
			return 0;
		}
	}
}

sub cons_position_is_in_high_quality_segment_of_read
{
	my ( $read, $cons_pos ) = @_;
	if ( $read->qual_clip_start == -1 && $read->qual_clip_stop == -1 )
	{
		return 0;
	}
	else
	{
		if ( $read->qual_clip_start <= $read->get_child_position_from_parent_position($cons_pos )
			&& $read->get_child_position_from_parent_position($cons_pos )<= $read->qual_clip_stop )
		{
			return 1;
		}
		else
		{
			return 0;
		}
	}
}


sub recalculate_consensus_qualities_and_change
{
	my ( $rContig, $nConsPosLeft, $nConsPosRight, $rReadArray ) = @_;

	my @aNewQualities;

	recalculate_consensus_qualities_no_change( $nConsPosLeft, $nConsPosRight, \@aNewQualities, $rReadArray, $rContig );

	my @qualities = @{$rContig->padded_base_quality};
	for ( my $nConsPos = $nConsPosLeft ; $nConsPos <= $nConsPosRight ; $nConsPos++ )
	{
		my $qualNew = $aNewQualities[ $nConsPos - $nConsPosLeft ];

		assert( QUALITY_MIN <= $qualNew );
		assert( $qualNew <= QUALITY_MAX );
		$qualities[$nConsPos - 1 ] = $qualNew;
	}
	$rContig->padded_base_quality(\@qualities);
}

sub check_contiguous
{

	my ( $reads, $start_position, $end_position ) = @_;

	my $pos;
	my @aCoveredByReadInRegion;
	for ( $pos = $start_position ; $pos <= $end_position ; $pos++ )
	{
		$aCoveredByReadInRegion[ $pos++ ] = 0;
	}

	foreach  my $read (@{$reads})
	{

		# I tried to decide whether to allow a tear when there might
		# only be unaligned sequence holding one of the resulting
		# contigs together.  I decided 'yes', because there would be
		# unaligned at the ends of each read, hence at the end of each
		# resulting contig.  That is normal, and should be allowed.

		my $read_start_position;
		my $read_end_position;

		assert(
			calculate_intersection( $start_position, $end_position, $read->align_start, $read->align_stop, \$read_start_position, \$read_end_position ) );

		for ( $pos = $read_start_position ; $pos <= $read_end_position ; $pos++ )
		{
			$aCoveredByReadInRegion[$pos] = 1;
		}
	}

	for ( $pos = $start_position ; $pos <= $end_position ; $pos++ )
	{
		if ( !$aCoveredByReadInRegion[$pos] )
		{

			# hole in the contig

			printf( "Contig has position %d that would not be covered by any read.\n", $pos );
			if ( wantarray() )
			{
				return ( 0, $pos );
			}
			else
			{
				return 0;
			}
		}
	}

	# if reached here, there is no hole in the part of the new contig in
	# the torn region
	if ( wantarray() )
	{
		return ( 1, $pos );
	}
	else
	{
		return 1;
	}
}

sub calculate_intersection
{
	my $bOK;
	my ( $nALeft, $nARight, $nBLeft, $nBRight, $nIntersectLeft, $nIntersectRight ) = @_;

	$$nIntersectLeft  = max( $nALeft,  $nBLeft );
	$$nIntersectRight = min( $nARight, $nBRight );

	if ( $$nIntersectLeft <= $$nIntersectRight )
	{
		$bOK = 1;
	}
	else
	{
		$bOK = 0;
	}
	return ($bOK);
}

sub intervals_intersect
{
	my ( $nL1, $nR1, $nL2, $nR2 ) = @_;
	return ( ( $nL1 <= $nR2 ) && ( $nL2 <= $nR1 ) );
}

sub recalculate_consensus_qualities_no_change
{
	my ( $nRegionConsPosLeft, $nRegionConsPosRight, $rNewQualitiesArray, $rReadArray, $rContig ) = @_;

	$#{$rNewQualitiesArray} = -1;
	my @aReadBoundaries;

	my $nGuessOfNumberOfBoundariesNeeded = ( ( $nRegionConsPosRight - $nRegionConsPosLeft ) * ( $#$rReadArray + 1 ) ) / $rContig->length;

	# the division site will be .5 base to the left of the number
	# put into this array

	push( @aReadBoundaries, $nRegionConsPosLeft );
	push( @aReadBoundaries, ( $nRegionConsPosRight + 1 ) );

	my $nRead;
	for ( $nRead = 0 ; $nRead < @$rReadArray ; $nRead++ )
	{
		my $rRead = $$rReadArray[$nRead];		

		my $nIntersectLeft;
		my $nIntersectRight;

		if (
			!calculate_intersection(
				$nRegionConsPosLeft, $nRegionConsPosRight,
				$rRead->get_parent_position_from_child_position( $rRead->align_clip_start ),
				$rRead->get_parent_position_from_child_position( $rRead->align_clip_stop ),
				\$nIntersectLeft, \$nIntersectRight
			)
		  )
		{
			next;
		}

		# if reached here, this read (partly) overlaps the region

		push( @aReadBoundaries, $nIntersectLeft );
		push( @aReadBoundaries, ( $nIntersectRight + 1 ) );
	}

	# now sort this.
	@aReadBoundaries = sort { $a <=> $b } @aReadBoundaries;

	# take the time and eliminate duplicates

	my $nSegment;

	for ( $nSegment = ( @aReadBoundaries - 1 ) ; $nSegment >= 1 ; $nSegment-- )
	{

		if ( $aReadBoundaries[$nSegment] == ( $aReadBoundaries[ $nSegment - 1 ] ) )
		{
			splice( @aReadBoundaries, $nSegment, 1 );
		}
	}

	my @aConsensusSegments;

	#go ahead and set the size of aConsensusSegments here so that perl doesn't have to keep resizing it
	for ( $nSegment = 0 ; $nSegment < ( @aReadBoundaries - 1 ) ; $nSegment++ )
	{

		my $nConsPosLeft  = $aReadBoundaries[$nSegment];
		my $nConsPosRight = $aReadBoundaries[ $nSegment + 1 ] - 1;

		push(
			@aConsensusSegments,
			{
				cons_pos_left  => $nConsPosLeft,
				cons_pos_right => $nConsPosRight
			}
		);
	}

	# go through all the reads again, this time putting each one into
	# each appropriate consensusSegment

	my @reads = map { {
	                                        name => $_->name, 
			                        align_clip_start => $_->align_clip_start, 
						align_clip_stop => $_->align_clip_stop, 
						padded_base_string => $_->padded_base_string,
						padded_base_quality => $_->padded_base_quality,
						position => $_->position,
						template => $_->template,
						qual_clip_start => $_->qual_clip_start,
						qual_clip_stop => $_->qual_clip_stop,
						complemented => $_->complemented,
						chem => $_->chemistry,
                        #chem => $_->chem,
						chromat_file => $_->chromat_file
						
					   } } @{$rReadArray};

	my $contig_padded_base_string = $rContig->padded_base_string;	

	for ( $nRead = 0 ; $nRead < @reads ; $nRead++ )
	{

		my $rRead = $reads[$nRead];

		my $nIntersectLeft;
		my $nIntersectRight;

		if (
			!calculate_intersection(
				$nRegionConsPosLeft, $nRegionConsPosRight,
				get_parent_position_from_child_position( $rRead, $rRead->{align_clip_start} ),
				get_parent_position_from_child_position( $rRead, $rRead->{align_clip_stop} ),
				\$nIntersectLeft, \$nIntersectRight
			)
		  )
		{
			next;
		}

		my $nFirstConsSeg = _nFindIndexOfMatchOrPredecessor( \@aReadBoundaries, $nIntersectLeft );

		# -------- ------- --------- ---------- consensus segments
		#    --------------------   read

		# ^ nFirstConsSeg

		my $nConsSeg;
		if ( !defined($nFirstConsSeg) )
		{
			$nConsSeg = 0;
		}
		else
		{
			$nConsSeg = $nFirstConsSeg;
		}

		my $bContinue = 1;

		for ( ; $nConsSeg < @aConsensusSegments && $bContinue ; $nConsSeg++ )
		{
			my $rConsSeg = $aConsensusSegments[$nConsSeg];

			if ( intervals_intersect( $nIntersectLeft, $nIntersectRight, $rConsSeg->{cons_pos_left}, $rConsSeg->{cons_pos_right} ) )
			{

				push( @{ $rConsSeg->{raReads} }, $rRead );
			}
			else
			{
				if ( $nIntersectRight < $rConsSeg->{cons_pos_left} )
				{

					# ------- ------- -------
					#   -----------   ^this consensusSegment
					#    read
					# So continuing to examine more consensusSegment 's will
					# not help since then will continue to be off the read
					# to the right.

					$bContinue = 0;
				}
			}
		}
	}

	my @aAgreeingReads;

	# put up here just so this array need not be created with each
	# consensus position

	my @aTemplates;

	# now work way through the consensusSegments, and within each of those,
	# through the consensus positions.  What is the highest quality read
	# at one base within a consensusSegment may be different than the highest
	# quality read at another base within the same consensusSegment.

	for ( $nSegment = 0 ; $nSegment < @aConsensusSegments ; $nSegment++ )
	{

		my $rConsSeg = $aConsensusSegments[$nSegment];

		for ( my $nConsPos = $rConsSeg->{cons_pos_left} ; $nConsPos <= $rConsSeg->{cons_pos_right} ; $nConsPos++ )
		{

			# we want to find a window about the consensus which includes
			# 5 non-pads.  I guess I will allow the consensus to be a pad
			# and just look for 2 non-pads on each side of the consensus base.

			# putting default values in case we are close to the beginning
			# or end of the sequence.  (I believe Purify found this problem
			# Jan 2002).
			my $nWindowLeft  = $nConsPos - 2;
			my $nWindowRight = $nConsPos + 2;

			my $nNonpads = 0;
			my $nConsPos3;
			for ( $nConsPos3 = $nConsPos - 1 ; $nConsPos3 >= 0 ; $nConsPos3-- )
			{
				if ( !( substr( $contig_padded_base_string, $nConsPos3 - 1, 1 ) eq '*' ) )    #if base at position $nConsPos3 is a pad then...
				{
					$nNonpads++;
					if ( $nNonpads == 2 )
					{
						$nWindowLeft = $nConsPos3;
						last;
					}
				}
			}

			$nNonpads = 0;
			my $contig_length = $rContig->length;
			for ( $nConsPos3 = $nConsPos + 1 ; $nConsPos3 < $contig_length ; $nConsPos3++ )
			{
				if ( !( substr( $contig_padded_base_string, $nConsPos3 - 1, 1 ) eq '*' ) )    #if base at position $nConsPos3 is a pad then...
				{
					$nNonpads++;
					if ( $nNonpads == 2 )
					{
						$nWindowRight = $nConsPos3;
						last;
					}
				}
			}

			# i need to fix this comment as soon as I figure out what it means, jks
			# consider reads that agrentGetFragFromConsPose with the
			# consensus in this window

			my $rBestAgreeingRead     = undef;
			my $nQualBestAgreeingRead = 0;
			@aTemplates = ();
			my $nNumberOfWholeCloneReads = 0;

			# we're looping the agreeing reads array, in consed, there are four types of read templates that we are looking for,
			# based on whether they belong to the top or bottom strand, and are formed by a dye primer or a dye terminator
			# these two variables can form four distinct kinds of "AgreeingReads" (i.e. bottom strand dye terminator, top strand dye primer...)
			# I don't think this loop is necessary, doing it the perl way instead
			for ( my $nGroup = 0 ; $nGroup < 4 ; $nGroup++ )
			{
				$aAgreeingReads[$nGroup] = [];
			}

			#$#aAgreeingReads=0;

			for ( my $nRead = 0 ; $nRead < @{ $rConsSeg->{raReads} } ; $nRead++ )
			{
				my $rRead = $rConsSeg->{raReads}[$nRead];

				# check first that the window is completely within the
				# aligned portion of the read

				if (
					!(
						get_parent_position_from_child_position( $rRead, $rRead->{align_clip_start} ) <= $nWindowLeft
						&& $nWindowRight <= get_parent_position_from_child_position( $rRead, $rRead->{align_clip_stop} )
					)
				  )
				{
					next;
				}

				my $bAgreeSoFar = 1;

				for ( my $nConsPos2 = $nWindowLeft ; $nConsPos2 <= $nWindowRight && $bAgreeSoFar ; $nConsPos2++ )
				{
					if (
						lc( substr( $rRead->{padded_base_string}, ( get_child_position_from_parent_position( $rRead, $nConsPos2 - 1 ) ), 1 ) ) ne
						lc( substr( $contig_padded_base_string, ( $nConsPos2 - 1 ), 1 ) ) )
					{
						$bAgreeSoFar = 0;
					}
				}

				if ($bAgreeSoFar)
				{
					my $nStrandAndChemistry = _nGetStrandAndChemistry($rRead);
					if ( !defined( $aAgreeingReads[$nStrandAndChemistry] ) )
					{
						$aAgreeingReads[$nStrandAndChemistry] = [];
					}
					push( @{ $aAgreeingReads[$nStrandAndChemistry] }, $rRead );

					my $nQual = $rRead->{padded_base_quality}->[ get_child_position_from_parent_position( $rRead, $nConsPos - 1 ) ];

					if ( !defined($rBestAgreeingRead)
						|| $nQual > $nQualBestAgreeingRead )
					{

						$nQualBestAgreeingRead = $nQual;
						$rBestAgreeingRead     = $rRead;
					}

					if ( _bIsWholeCloneTemplateNotSubcloneTemplate($rRead) )
					{
						$nNumberOfWholeCloneReads++;
					}
					else
					{
						my $bMatchFound = 0;
						for ( my $nIndex = 0 ; $nIndex < @aTemplates ; $nIndex++ )
						{
							if ( $aTemplates[$nIndex] eq $rRead->{template} )
							{
								$bMatchFound = 1;
							}
						}
						push( @aTemplates, $rRead->{template} )
						  if !$bMatchFound;
					}
				}
			}

			# we've already found the highest quality *single* read (not part
			# of a pair)

			my $rBestRead1 = $rBestAgreeingRead;
			my $rBestRead2;
			my $nQualBest = $nQualBestAgreeingRead;

			# now look over pairs of reads to see if we can
			# find a pair of reads that does better

			# pick reads from each of the 1st 3 groups
			for ( my $nStrandChem1 = 0 ; $nStrandChem1 <= 2 ; $nStrandChem1++ )
			{

				# pick all the reads in the group
				for ( my $nRead1 = 0 ; defined( $aAgreeingReads[$nStrandChem1] ) && $nRead1 < @{ $aAgreeingReads[$nStrandChem1] } ; $nRead1++ )
				{
					my $rReadStrandChem1 = $aAgreeingReads[$nStrandChem1][$nRead1];

					# pick reads from each of the remaining groups
					for ( my $nStrandChem2 = $nStrandChem1 + 1 ; $nStrandChem2 <= 3 ; $nStrandChem2++ )
					{

						# pick all reads in the group
						for ( my $nRead2 = 0 ; $nRead2 < ( $#{ $aAgreeingReads[$nStrandChem2] } + 1 ) ; $nRead2++ )
						{
							my $rReadStrandChem2 = $aAgreeingReads[$nStrandChem2][$nRead2];

							next
							  if ( _bAreTheseReadsFromTheSameTemplate( $rReadStrandChem1, $rReadStrandChem2 ) );

							# if reach here, the reads are from different
							# templates

							my $nQual =
							  $rReadStrandChem1->{padded_base_quality}->[ get_child_position_from_parent_position( $rReadStrandChem1, $nConsPos - 1 ) ] +
							  $rReadStrandChem2->{padded_base_quality}->[ get_child_position_from_parent_position( $rReadStrandChem2, $nConsPos - 1 ) ];

							if ( $nQual > $nQualBest )
							{
								$rBestRead1 = $rReadStrandChem1;
								$rBestRead2 = $rReadStrandChem2;
								$nQualBest  = $nQual;

							}
						}    #  for( $nRead2 = 0;
					}    #  for( $nStrandChem2 ...
				}    #  for( $nRead1 = 0;
			}    #  for( $nStrandChem1 ...

			# when reached here, we have tried each individual read
			# and we have tried all pairs of reads from different
			# chem/strand groups.  The pair pBestRead1, pBestRead2
			# represents the best pair or best single read ( if pBestRead2
			# is null).

			# check if we quality for a +5 boost.  If there is a pair
			# of reads that gives the best quality, then this pair already
			# has 2 templates, so there must be 3 templates or more to
			# give a +5 boost.  If there is a single read, then there must
			# be 2 or more templates to give a +5 boost.

			my $nNumberOfTemplates = $nNumberOfWholeCloneReads + @aTemplates;

			if ($rBestRead2)
			{
				if ( $nNumberOfTemplates >= 3 )
				{
					$nQualBest += 5;
				}
			}
			else
			{
				if ( $nNumberOfTemplates >= 2 )
				{
					$nQualBest += 5;
				}
			}

			# there is an exception.  If there is only a single template,
			# and all the reads are in the low quality segment, then
			# set the consensus to quality 0.

			# The reason we don't check for == 0 is that in this case the
			# quality will already be set to zero.  The reason we don't need
			# to check the case of nNumberOfTemplates > 1 is that in this case
			# there are > 1 aligned templates perfectly matching a 5 base window,
			# so there are clearly > 1 aligned templates

			if ( $nNumberOfTemplates == 1 )
			{

				# are all of the reads in the low quality segment?

				my $bAllInLowQualitySegment = 1;

				# we are just interested in aligned reads--not whether
				# they agree with the consensus in a 5 base window
				for ( my $nRead = 0 ; $nRead < @{ $rConsSeg->{raReads} } && $bAllInLowQualitySegment ; $nRead++ )
				{
					my $rRead = $rConsSeg->{raReads}[$nRead];

					if ( _cons_position_is_in_high_quality_segment_of_read( $rRead, $nConsPos ) )
					{
						$bAllInLowQualitySegment = 0;
					}
				}

				if ($bAllInLowQualitySegment)
				{

					# now check if there is only one template that is aligned
					# When we checked before nNumberOfTemplates == 1, we
					# were counting the number of templates with 5 base
					# window agreeing with the consensus.  Now let's do
					# a more careful check that there is only one that is
					# aligned:

					my @aAlignedTemplates;
					my $bOnlyOneAlignedTemplate = 1;

					for ( my $nRead = 0 ; $nRead < ( $#{ $rConsSeg->{raReads} } + 1 ) && $bOnlyOneAlignedTemplate ; $nRead++ )
					{
						my $rRead       = $rConsSeg->{raReads}[$nRead];
						my $bMatchFound = 0;
						if ( _bIsWholeCloneTemplateNotSubcloneTemplate($rRead) )
						{
							$bOnlyOneAlignedTemplate = 0;
						}
						else
						{
							for ( my $nIndex = 0 ; $nIndex < @aAlignedTemplates ; $nIndex++ )
							{
								if ( $aAlignedTemplates[$nIndex] eq $rRead->{template} )
								{
									$bMatchFound = 1;
								}
							}
						}
						if ( !$bMatchFound )
						{
							push( @aAlignedTemplates, $rRead->{template} );
							if ( @aAlignedTemplates > 1 )
							{
								$bOnlyOneAlignedTemplate = 0;
							}

						}
					}

					if ($bOnlyOneAlignedTemplate)
					{
						$nQualBest = 0;
					}
				}    #  if ( $bAllInLowQualitySegment ) {
			}    #  if ( $nNumberOfTemplates == 1 ) {

			# now we have found the answer:  nQualBest.  We need to put the
			# answer in the Sequence

			if ( $nQualBest > 90 )
			{
				$nQualBest = 90;
			}			
			push( @{$rNewQualitiesArray}, $nQualBest );

		}    # nConsPos
	}    # for( int nSegment = 0; nSegment < aConsensusSegments.length();

}

sub get_child_position_from_parent_position
{
	my ($self, $parent_pos) = @_;
	return ( $parent_pos + ( 1 - $self->{position} ) );    #removed -1 at the end
}

sub get_parent_position_from_child_position
{
	my ($self, $child_pos) = @_;
	return ( ( $self->{position} - 1 ) + $child_pos );    #removed +1 at the end
}

sub recalculate_consensus_qualities_no_change_old
{
	my ( $nRegionConsPosLeft, $nRegionConsPosRight, $rNewQualitiesArray, $rReadArray, $rContig ) = @_;

	$#{$rNewQualitiesArray} = -1;
	my @aReadBoundaries;

	my $nGuessOfNumberOfBoundariesNeeded = ( ( $nRegionConsPosRight - $nRegionConsPosLeft ) * ( $#$rReadArray + 1 ) ) / $rContig->length;

	# the division site will be .5 base to the left of the number
	# put into this array

	push( @aReadBoundaries, $nRegionConsPosLeft );
	push( @aReadBoundaries, ( $nRegionConsPosRight + 1 ) );

	my $nRead;
	for ( $nRead = 0 ; $nRead < @$rReadArray ; $nRead++ )
	{
		my $rRead = $$rReadArray[$nRead];		

		my $nIntersectLeft;
		my $nIntersectRight;

		if (
			!calculate_intersection(
				$nRegionConsPosLeft, $nRegionConsPosRight,
				$rRead->get_parent_position_from_child_position( $rRead->align_clip_start ),
				$rRead->get_parent_position_from_child_position( $rRead->align_clip_stop ),
				\$nIntersectLeft, \$nIntersectRight
			)
		  )
		{
			next;
		}

		# if reached here, this read (partly) overlaps the region

		push( @aReadBoundaries, $nIntersectLeft );
		push( @aReadBoundaries, ( $nIntersectRight + 1 ) );
	}

	# now sort this.
	@aReadBoundaries = sort { $a <=> $b } @aReadBoundaries;

	# take the time and eliminate duplicates

	my $nSegment;

	for ( $nSegment = ( @aReadBoundaries - 1 ) ; $nSegment >= 1 ; $nSegment-- )
	{

		if ( $aReadBoundaries[$nSegment] == ( $aReadBoundaries[ $nSegment - 1 ] ) )
		{
			splice( @aReadBoundaries, $nSegment, 1 );
		}
	}

	my @aConsensusSegments;

	#go ahead and set the size of aConsensusSegments here so that perl doesn't have to keep resizing it
	for ( $nSegment = 0 ; $nSegment < ( @aReadBoundaries - 1 ) ; $nSegment++ )
	{

		my $nConsPosLeft  = $aReadBoundaries[$nSegment];
		my $nConsPosRight = $aReadBoundaries[ $nSegment + 1 ] - 1;

		push(
			@aConsensusSegments,
			{
				cons_pos_left  => $nConsPosLeft,
				cons_pos_right => $nConsPosRight
			}
		);
	}

	# go through all the reads again, this time putting each one into
	# each appropriate consensusSegment

	for ( $nRead = 0 ; $nRead < @$rReadArray ; $nRead++ )
	{

		my $rRead = $rReadArray->[$nRead];

		my $nIntersectLeft;
		my $nIntersectRight;

		if (
			!calculate_intersection(
				$nRegionConsPosLeft, $nRegionConsPosRight,
				$rRead->get_parent_position_from_child_position( $rRead->align_clip_start ),
				$rRead->get_parent_position_from_child_position( $rRead->align_clip_stop ),
				\$nIntersectLeft, \$nIntersectRight
			)
		  )
		{
			next;
		}

		my $nFirstConsSeg = _nFindIndexOfMatchOrPredecessor( \@aReadBoundaries, $nIntersectLeft );

		# -------- ------- --------- ---------- consensus segments
		#    --------------------   read

		# ^ nFirstConsSeg

		my $nConsSeg;
		if ( !defined($nFirstConsSeg) )
		{
			$nConsSeg = 0;
		}
		else
		{
			$nConsSeg = $nFirstConsSeg;
		}

		my $bContinue = 1;

		for ( ; $nConsSeg < @aConsensusSegments && $bContinue ; $nConsSeg++ )
		{
			my $rConsSeg = $aConsensusSegments[$nConsSeg];

			if ( intervals_intersect( $nIntersectLeft, $nIntersectRight, $rConsSeg->{cons_pos_left}, $rConsSeg->{cons_pos_right} ) )
			{

				push( @{ $rConsSeg->{raReads} }, $rRead );
			}
			else
			{
				if ( $nIntersectRight < $rConsSeg->{cons_pos_left} )
				{

					# ------- ------- -------
					#   -----------   ^this consensusSegment
					#    read
					# So continuing to examine more consensusSegment 's will
					# not help since then will continue to be off the read
					# to the right.

					$bContinue = 0;
				}
			}
		}
	}

	my @aAgreeingReads;

	# put up here just so this array need not be created with each
	# consensus position

	my @aTemplates;

	# now work way through the consensusSegments, and within each of those,
	# through the consensus positions.  What is the highest quality read
	# at one base within a consensusSegment may be different than the highest
	# quality read at another base within the same consensusSegment.

	for ( $nSegment = 0 ; $nSegment < @aConsensusSegments ; $nSegment++ )
	{

		my $rConsSeg = $aConsensusSegments[$nSegment];

		for ( my $nConsPos = $rConsSeg->{cons_pos_left} ; $nConsPos <= $rConsSeg->{cons_pos_right} ; $nConsPos++ )
		{

			# we want to find a window about the consensus which includes
			# 5 non-pads.  I guess I will allow the consensus to be a pad
			# and just look for 2 non-pads on each side of the consensus base.

			# putting default values in case we are close to the beginning
			# or end of the sequence.  (I believe Purify found this problem
			# Jan 2002).
			my $nWindowLeft  = $nConsPos - 2;
			my $nWindowRight = $nConsPos + 2;

			my $nNonpads = 0;
			my $nConsPos3;
			for ( $nConsPos3 = $nConsPos - 1 ; $nConsPos3 >= 0 ; $nConsPos3-- )
			{
				if ( !( substr( $rContig->padded_base_string, $nConsPos3 - 1, 1 ) eq '*' ) )    #if base at position $nConsPos3 is a pad then...
				{
					$nNonpads++;
					if ( $nNonpads == 2 )
					{
						$nWindowLeft = $nConsPos3;
						last;
					}
				}
			}

			$nNonpads = 0;
			for ( $nConsPos3 = $nConsPos + 1 ; $nConsPos3 < $rContig->length ; $nConsPos3++ )
			{
				if ( !( substr( $rContig->padded_base_string, $nConsPos3 - 1, 1 ) eq '*' ) )    #if base at position $nConsPos3 is a pad then...
				{
					$nNonpads++;
					if ( $nNonpads == 2 )
					{
						$nWindowRight = $nConsPos3;
						last;
					}
				}
			}

			# i need to fix this comment as soon as I figure out what it means, jks
			# consider reads that agrentGetFragFromConsPose with the
			# consensus in this window

			my $rBestAgreeingRead     = undef;
			my $nQualBestAgreeingRead = 0;
			@aTemplates = ();
			my $nNumberOfWholeCloneReads = 0;

			# we're looping the agreeing reads array, in consed, there are four types of read templates that we are looking for,
			# based on whether they belong to the top or bottom strand, and are formed by a dye primer or a dye terminator
			# these two variables can form four distinct kinds of "AgreeingReads" (i.e. bottom strand dye terminator, top strand dye primer...)
			# I don't think this loop is necessary, doing it the perl way instead
			for ( my $nGroup = 0 ; $nGroup < 4 ; $nGroup++ )
			{
				$aAgreeingReads[$nGroup] = [];
			}

			#$#aAgreeingReads=0;

			for ( my $nRead = 0 ; $nRead < @{ $rConsSeg->{raReads} } ; $nRead++ )
			{
				my $rRead = $rConsSeg->{raReads}[$nRead];

				# check first that the window is completely within the
				# aligned portion of the read

				if (
					!(
						$rRead->get_parent_position_from_child_position( $rRead->align_clip_start ) <= $nWindowLeft
						&& $nWindowRight <= $rRead->get_parent_position_from_child_position( $rRead->align_clip_stop )
					)
				  )
				{
					next;
				}

				my $bAgreeSoFar = 1;

				for ( my $nConsPos2 = $nWindowLeft ; $nConsPos2 <= $nWindowRight && $bAgreeSoFar ; $nConsPos2++ )
				{
					if (
						lc( substr( $rRead->padded_base_string, ( $rRead->get_child_position_from_parent_position( $nConsPos2 - 1 ) ), 1 ) ) ne
						lc( substr( $rContig->padded_base_string, ( $nConsPos2 - 1 ), 1 ) ) )
					{
						$bAgreeSoFar = 0;
					}
				}

				if ($bAgreeSoFar)
				{
					my $nStrandAndChemistry = _nGetStrandAndChemistry($rRead);
					if ( !defined( $aAgreeingReads[$nStrandAndChemistry] ) )
					{
						$aAgreeingReads[$nStrandAndChemistry] = [];
					}
					push( @{ $aAgreeingReads[$nStrandAndChemistry] }, $rRead );

					my $nQual = $rRead->padded_base_quality->[ $rRead->get_child_position_from_parent_position( $nConsPos - 1 ) ];

					if ( !defined($rBestAgreeingRead)
						|| $nQual > $nQualBestAgreeingRead )
					{

						$nQualBestAgreeingRead = $nQual;
						$rBestAgreeingRead     = $rRead;
					}

					if ( _bIsWholeCloneTemplateNotSubcloneTemplate($rRead) )
					{
						$nNumberOfWholeCloneReads++;
					}
					else
					{
						my $bMatchFound = 0;
						for ( my $nIndex = 0 ; $nIndex < @aTemplates ; $nIndex++ )
						{
							if ( $aTemplates[$nIndex] eq $rRead->template )
							{
								$bMatchFound = 1;
							}
						}
						push( @aTemplates, $rRead->template )
						  if !$bMatchFound;
					}
				}
			}

			# we've already found the highest quality *single* read (not part
			# of a pair)

			my $rBestRead1 = $rBestAgreeingRead;
			my $rBestRead2;
			my $nQualBest = $nQualBestAgreeingRead;

			# now look over pairs of reads to see if we can
			# find a pair of reads that does better

			# pick reads from each of the 1st 3 groups
			for ( my $nStrandChem1 = 0 ; $nStrandChem1 <= 2 ; $nStrandChem1++ )
			{

				# pick all the reads in the group
				for ( my $nRead1 = 0 ; defined( $aAgreeingReads[$nStrandChem1] ) && $nRead1 < @{ $aAgreeingReads[$nStrandChem1] } ; $nRead1++ )
				{
					my $rReadStrandChem1 = $aAgreeingReads[$nStrandChem1][$nRead1];

					# pick reads from each of the remaining groups
					for ( my $nStrandChem2 = $nStrandChem1 + 1 ; $nStrandChem2 <= 3 ; $nStrandChem2++ )
					{

						# pick all reads in the group
						for ( my $nRead2 = 0 ; $nRead2 < ( $#{ $aAgreeingReads[$nStrandChem2] } + 1 ) ; $nRead2++ )
						{
							my $rReadStrandChem2 = $aAgreeingReads[$nStrandChem2][$nRead2];

							next
							  if ( _bAreTheseReadsFromTheSameTemplate( $rReadStrandChem1, $rReadStrandChem2 ) );

							# if reach here, the reads are from different
							# templates

							my $nQual =
							  $rReadStrandChem1->padded_base_quality->[ $rReadStrandChem1->get_child_position_from_parent_position( $nConsPos - 1 ) ] +
							  $rReadStrandChem2->padded_base_quality->[ $rReadStrandChem2->get_child_position_from_parent_position( $nConsPos - 1 ) ];

							if ( $nQual > $nQualBest )
							{
								$rBestRead1 = $rReadStrandChem1;
								$rBestRead2 = $rReadStrandChem2;
								$nQualBest  = $nQual;

							}
						}    #  for( $nRead2 = 0;
					}    #  for( $nStrandChem2 ...
				}    #  for( $nRead1 = 0;
			}    #  for( $nStrandChem1 ...

			# when reached here, we have tried each individual read
			# and we have tried all pairs of reads from different
			# chem/strand groups.  The pair pBestRead1, pBestRead2
			# represents the best pair or best single read ( if pBestRead2
			# is null).

			# check if we quality for a +5 boost.  If there is a pair
			# of reads that gives the best quality, then this pair already
			# has 2 templates, so there must be 3 templates or more to
			# give a +5 boost.  If there is a single read, then there must
			# be 2 or more templates to give a +5 boost.

			my $nNumberOfTemplates = $nNumberOfWholeCloneReads + @aTemplates;

			if ($rBestRead2)
			{
				if ( $nNumberOfTemplates >= 3 )
				{
					$nQualBest += 5;
				}
			}
			else
			{
				if ( $nNumberOfTemplates >= 2 )
				{
					$nQualBest += 5;
				}
			}

			# there is an exception.  If there is only a single template,
			# and all the reads are in the low quality segment, then
			# set the consensus to quality 0.

			# The reason we don't check for == 0 is that in this case the
			# quality will already be set to zero.  The reason we don't need
			# to check the case of nNumberOfTemplates > 1 is that in this case
			# there are > 1 aligned templates perfectly matching a 5 base window,
			# so there are clearly > 1 aligned templates

			if ( $nNumberOfTemplates == 1 )
			{

				# are all of the reads in the low quality segment?

				my $bAllInLowQualitySegment = 1;

				# we are just interested in aligned reads--not whether
				# they agree with the consensus in a 5 base window
				for ( my $nRead = 0 ; $nRead < @{ $rConsSeg->{raReads} } && $bAllInLowQualitySegment ; $nRead++ )
				{
					my $rRead = $rConsSeg->{raReads}[$nRead];

					if ( cons_position_is_in_high_quality_segment_of_read( $rRead, $nConsPos ) )
					{
						$bAllInLowQualitySegment = 0;
					}
				}

				if ($bAllInLowQualitySegment)
				{

					# now check if there is only one template that is aligned
					# When we checked before nNumberOfTemplates == 1, we
					# were counting the number of templates with 5 base
					# window agreeing with the consensus.  Now let's do
					# a more careful check that there is only one that is
					# aligned:

					my @aAlignedTemplates;
					my $bOnlyOneAlignedTemplate = 1;

					for ( my $nRead = 0 ; $nRead < ( $#{ $rConsSeg->{raReads} } + 1 ) && $bOnlyOneAlignedTemplate ; $nRead++ )
					{
						my $rRead       = $rConsSeg->{raReads}[$nRead];
						my $bMatchFound = 0;
						if ( _bIsWholeCloneTemplateNotSubcloneTemplate($rRead) )
						{
							$bOnlyOneAlignedTemplate = 0;
						}
						else
						{
							for ( my $nIndex = 0 ; $nIndex < @aAlignedTemplates ; $nIndex++ )
							{
								if ( $aAlignedTemplates[$nIndex] eq $rRead->template )
								{
									$bMatchFound = 1;
								}
							}
						}
						if ( !$bMatchFound )
						{
							push( @aAlignedTemplates, $rRead->template );
							if ( @aAlignedTemplates > 1 )
							{
								$bOnlyOneAlignedTemplate = 0;
							}

						}
					}

					if ($bOnlyOneAlignedTemplate)
					{
						$nQualBest = 0;
					}
				}    #  if ( $bAllInLowQualitySegment ) {
			}    #  if ( $nNumberOfTemplates == 1 ) {

			# now we have found the answer:  nQualBest.  We need to put the
			# answer in the Sequence

			if ( $nQualBest > 90 )
			{
				$nQualBest = 90;
			}			
			push( @{$rNewQualitiesArray}, $nQualBest );

		}    # nConsPos
	}    # for( int nSegment = 0; nSegment < aConsensusSegments.length();

}

1;
