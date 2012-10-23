#The Contig Merging and Splitting Utility Library has several classes.
#The first class is a general
package Utility;
use Exporter;
use Compress::Zlib;
@ISA    = qw(Exporter);
@EXPORT = qw(
@aUnPaddedToPadded 
@aPaddedToUnPadded 
max
min
bIntersect
bIntervalsIntersect
nFindIndexOfMatchOrPredecessor

CreatePaddedAndUnPaddedArrays
CreatePaddedConsQualityArray
CreateUnPaddedConsQualityArray
SetConsensusToProperCapitals
recalculateConsensusQualitiesNoChange
recalculateConsensusQualitiesAndChange 
bCheckContiguous

CreatePaddedQualityArray	
SetSequenceToProperCapitals
nNormalQualityFrom9899Quality
nGetReadIndexFromConsPos
nGetConsPosFromReadIndex
nGetAlignStart
nGetAlignEnd

_nGetStrandAndChemistry
_bIsDyePrimerNotDyeTerminator
_bComp
_bIsWholeCloneTemplateNotSubcloneTemplate
_bIsInHighQualitySegmentOfRead	

bGetBSArrayStructureOk	
nGetBaseSegmentIndexByPos
numerically

LoadPhdInfo
);

our $VERSION = 0.01;

use strict;

#use warnings;
use Carp;
use Carp::Assert;
use Cwd;
use IO::String;

our @aUnPaddedToPadded;
our @aPaddedToUnPadded;

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

sub LoadPhdInfo
{
	my ( $raRead, $nLeftPos, $nRightPos, $sGetPhdFrom ) = @_;
	#foreach my $tempread (@$raRead)
	#{
	#	print "$tempread->{name}\n";
	#}
	my %phdreads;     #hash containing the required phd files as the key, and it's
	                  #corresponding read, we check against this when loading fasta.qual and phd
	                  #files below.
	my %readreads;    #hash of reads intersecting the merge region

	if ( $sGetPhdFrom eq "DB" )    #get qual info from database
	{

		my @overlapping_read_names;

		#first we build our readread hash table
		for ( my $nPos = 0 ; $nPos < @{$raRead} ; $nPos++ )
		{
			my $rRead = ${$raRead}[$nPos];
			if ( bIntervalsIntersect( nGetAlignStart($rRead), nGetAlignEnd($rRead), $nLeftPos, $nRightPos ) )
			{
				my $rPhd = {};
				if ( !defined( $rPhd->{template} ) )
				{
					my $tempname = $rRead->{name};
					$tempname =~ s/(\..+)$//;    #chop off the extension and use it to define our template type
					$rPhd->{template} = $tempname;
				}
				$rRead->{phd} = $rPhd;
				$readreads{ $rRead->{name} } = $rRead;
				my ($phdversion) = $rRead->{description}{PHD_FILE} =~ /\.(\d+)$/;
				push( @overlapping_read_names, $rRead->{name} . "-$phdversion" );
			}
		}

		my @overlapping_reads = GSC::Sequence::Item->get( sequence_item_name => \@overlapping_read_names );

		my @read_ids = map { $_->seq_id } @overlapping_reads;
		my @phds = GSC::Sequence::PHD->get( seq_id => \@read_ids );

		#  #need to get quality values for each read
		my $rPhdReader = new GSC::IO::Assembly::Phd::Reader();
		my $sPHDFILE;
		my $rPhd;
		my $read;
		foreach $read (@overlapping_reads)
		{
			my ($read_name) = $read->sequence_item_name =~ /(.+)-\d+$/;
			my $rRead       = $readreads{$read_name};
			my $phd_content = $read->phd_content;
			my $phd_file = new IO::String($phd_content);

			my $rPhd = $rPhdReader->read($phd_file);

			if ( !defined( $rPhd->{template} ) )
			{
				my $tempname = $rRead->{name};
				$tempname =~ s/(\..+)$//;    #chop off the extension and use it to define our template type
				$rPhd->{template} = $tempname;
			}

			if ( _bComp($rRead) )
			{
				@{ $rPhd->{sequence}{quality} } = reverse @{ $rPhd->{sequence}{quality} };
			}

			$rRead->{phd} = $rPhd;
			CreatePaddedQualityArray($rRead);
			$phd_file->close;
		}
	}
	elsif ( -e "../phd_dir" )
	{

		#need to get quality values for each read
		my $sPWD       = cwd();
		my $sPHDDIR    = $sPWD . "/../phd_dir/";
		my $rPhdReader = new GSC::IO::Assembly::Phd::Reader();
		my $sPHDFILE;
		my $rPhd;

		for ( my $nPos = 0 ; $nPos < @{$raRead} ; $nPos++ )
		{
			my $rRead = ${$raRead}[$nPos];
			if ( bIntervalsIntersect( nGetAlignStart($rRead), nGetAlignEnd($rRead), $nLeftPos, $nRightPos ) )
			{
				$sPHDFILE = $sPHDDIR . $rRead->{description}->{PHD_FILE};
				my $PHDFILE;
				if ( -e $sPHDFILE )
				{
					$PHDFILE = new IO::File($sPHDFILE);    #
					                                       #open( PHDFILE, $sPHDFILE );
				}
				elsif ( -e "$sPHDFILE.gz" )
				{

					#	$PHDFILE = Compress::Zlib::gzopen("$sPHDFILE.gz", "ro");
					open( PHDFILE, "$sPHDFILE.gz" );
					my $sPhdContent = join( '', <PHDFILE> );
					my $PHDFILEDATA = Compress::Zlib::memGunzip($sPhdContent);
					$PHDFILE = new IO::String($PHDFILEDATA);
				}

				my $rPhd = $rPhdReader->read($PHDFILE);

				if ( !defined( $rPhd->{template} ) )
				{
					my $tempname = $rRead->{name};
					$tempname =~ s/(\..+)$//;    #chop off the extension and use it to define our template type
					$rPhd->{template} = $tempname;
				}

				if ( _bComp($rRead) )
				{
					@{ $rPhd->{sequence}{quality} } = reverse @{ $rPhd->{sequence}{quality} };
				}

				$rRead->{phd} = $rPhd;
				CreatePaddedQualityArray($rRead);
				close(PHDFILE);
			}
		}
	}
	elsif ( opendir THISDIR, "../Input" )
	{
		my @allFiles  = readdir THISDIR;
		my $bUsingTab = 0;
		my $sTabFile;
		foreach (@allFiles)
		{
			my ( $Name, $Ext ) = $_ =~ /^(.+\.)(.+)$/;
			if ( "tab" eq $Ext )
			{
				$bUsingTab = 1;
				$sTabFile  = $_;
				print "Using Tab File as read index.\n";
				last;
			}
		}
		if ($bUsingTab)    #we're not going to create hash since it's faster to do it on the fly
		{
			my $line;
			my $qualline;

			#first we build our readread hash table
			for ( my $nPos = 0 ; $nPos < @{$raRead} ; $nPos++ )
			{
				my $rRead = ${$raRead}[$nPos];
				if ( bIntervalsIntersect( nGetAlignStart($rRead), nGetAlignEnd($rRead), $nLeftPos, $nRightPos ) )
				{
					my $rPhd = {};
					if ( !defined( $rPhd->{template} ) )
					{
						my $tempname = $rRead->{name};
						$tempname =~ s/(\..+)$//;    #chop off the extension and use it to define our template type
						$rPhd->{template} = $tempname;
					}
					$rRead->{phd} = $rPhd;
					$readreads{ $rRead->{name} } = $rRead;
				}
			}
			$sTabFile = "../Input/" . $sTabFile;
			print "$sTabFile \n";
			open( TABFILE, "$sTabFile" ) || print "Failed to open $sTabFile\n";
			my $btabDontReadLine = 0;
			while ( $btabDontReadLine || ( $line = <TABFILE> ) )
			{
				$btabDontReadLine = 0;
				if ( $line =~ /\>/ )
				{
					my $sQualFileName = substr( $line, 2, ( length $line ) - 2 );
					chomp($sQualFileName);

					my $temp = substr( $line, 2, ( length($line) ) - 2 );
					chomp($temp);
					$sQualFileName = "../Input/" . $sQualFileName;
					unless ( open( QUALFILE, $sQualFileName ) )
					{
						print "Failed to open $sQualFileName\n";
						print "$line \n";
						print "$temp \n";
					}

				}
				else
				{
					next;
				}
				while ( $line = <TABFILE> )
				{
					if ( $line =~ /\>/ )
					{
						$btabDontReadLine = 1;
						last;
					}
					my @sTabTokens = split( / /, $line );
					my $sReadName = shift(@sTabTokens);
					chomp $sReadName;
					if ( exists $readreads{$sReadName} )
					{

						my $nQualFilePos = shift(@sTabTokens);
						chomp $nQualFilePos;
						my $rRead = $readreads{$sReadName};
						my @bq;
						seek QUALFILE, $nQualFilePos, 0;
						$qualline = <QUALFILE>;    #go to the next line
						while ( $qualline = <QUALFILE> )
						{

							if ( $qualline =~ /\>/ )
							{
								last;
							}
							chomp $qualline;
							push @bq, split( / /, $qualline );
						}
						$rRead->{phd}{sequence} = {};
						$rRead->{phd}{sequence}{quality} = \@bq;
						if ( _bComp($rRead) )
						{
							@{ $rRead->{phd}{sequence}{quality} } = reverse @{ $rRead->{phd}{sequence}{quality} };
						}
						CreatePaddedQualityArray($rRead);
					}
				}
				close QUALFILE;
			}
			return;

		}

		#if we reach here then we're doing it the slow way...

		#first we build our readread hash table
		for ( my $nPos = 0 ; $nPos < @{$raRead} ; $nPos++ )
		{
			my $rRead = ${$raRead}[$nPos];
			if ( bIntervalsIntersect( nGetAlignStart($rRead), nGetAlignEnd($rRead), $nLeftPos, $nRightPos ) )
			{
				my $rPhd = {};
				if ( !defined( $rPhd->{template} ) )
				{
					my $tempname = $rRead->{name};
					$tempname =~ s/(\..+)$//;    #chop off the extension and use it to define our template type
					$rPhd->{template} = $tempname;
				}
				$rRead->{phd} = $rPhd;
				$readreads{ $rRead->{name} } = $rRead;
			}
		}
		foreach (@allFiles)
		{
			my ( $Name, $Ext ) = $_ =~ /^(.+\.)(.+)$/;
			next unless "qual" eq $Ext;
			my $sQualFile = "../Input/" . $_;
			if ( open( QUALFILE, $sQualFile ) )
			{
				my $line;

				my $bDontReadLine = 0;    #set this if we don't need to read the line
				while ( $bDontReadLine || ( $line = <QUALFILE> ) )
				{
					$bDontReadLine = 0;
					next unless ( $line =~ /\>/ );
					my @tokens = split( / /, $line );
					my $readfilename = shift @tokens;
					chomp $readfilename;
					$readfilename = substr( $readfilename, 1, ( length($readfilename) - 1 ) );
					if ( exists $readreads{$readfilename} )
					{
						my $rRead = $readreads{$readfilename};
						my @bq;

						while ( $line = <QUALFILE> )
						{
							if ( $line =~ /\>/ )
							{
								$bDontReadLine = 1;

								last;
							}
							chomp $line;
							push @bq, split( / /, $line );
						}
						$rRead->{phd}{sequence} = {};
						$rRead->{phd}{sequence}{quality} = \@bq;
						if ( _bComp($rRead) )
						{
							@{ $rRead->{phd}{sequence}{quality} } = reverse @{ $rRead->{phd}{sequence}{quality} };
						}
						CreatePaddedQualityArray($rRead);
					}
				}
			}
		}
	}
	else
	{
		die
"Could not find a valid source of phd information.\n  The cmt needs to get phd information from either the database, an ../Input directory, or a ../phd_dir.\n";
	}

}

sub CreatePaddedAndUnPaddedArrays
{
	my ($rContig) = @_;
	my $nUnPadIndex = 0;
	for ( my $nPos = 0 ; $nPos < length( $rContig->{consensus} ) ; $nPos++ )
	{
		$aPaddedToUnPadded[$nPos]        = $nUnPadIndex;
		$aUnPaddedToPadded[$nUnPadIndex] = $nPos;
		if ( substr( $rContig->{consensus}, $nPos, 1 ) ne '*' )
		{
			$nUnPadIndex++;
		}
	}
}

sub nUnPaddedToPadded 
{
	my ($rContig, $nUnPaddedPosition) = @_;
	return $aUnPaddedToPadded[$nUnPaddedPosition];
}

sub nPaddedToUnPadded 
{
	my ($rContig, $nPaddedPosition) = @_;
	return $aPaddedToUnPadded[$nPaddedPosition];
}

sub CreatePaddedQualityArray
{
	my ($rRead)       = @_;
	my $rNewQualArray = [];
	my $nQualIndex    = 0;
	if ( !defined( $rRead->{phd} ) || !defined( $rRead->{phd}{sequence} ) || !defined( $rRead->{phd}{sequence}{quality} ) ) { return 1; }
	for ( my $nPos = 0 ; $nPos < length( $rRead->{sequence} ) ; $nPos++ )
	{
		if ( substr( $rRead->{sequence}, $nPos, 1 ) ne '*' )
		{
			$rNewQualArray->[$nPos] = $rRead->{phd}{sequence}{quality}[$nQualIndex];
			$nQualIndex++;
		}
		else
		{
			if ( defined( $rRead->{phd}{sequence}{quality}[ $nQualIndex - 1 ] ) && defined( $rRead->{phd}{sequence}{quality}[ $nQualIndex + 1 ] ) )
			{
				$rNewQualArray->[$nPos] =
				  int( ( $rRead->{phd}{sequence}{quality}[ $nQualIndex - 1 ] + $rRead->{phd}{sequence}{quality}[ $nQualIndex + 1 ] ) / 2 );
			}
		}
	}
	$rRead->{phd}{sequence}{quality} = $rNewQualArray;

	#yes, I know this is nasty, but it's the only quick fix I could think of
	#I'll create a separate data container later
}

sub CreatePaddedConsQualityArray
{
	my ($rContig)     = @_;
	my $rNewQualArray = [];
	my $nQualIndex    = 0;
	for ( my $nPos = 0 ; $nPos < length( $rContig->{consensus} ) ; $nPos++ )
	{
		if ( substr( $rContig->{consensus}, $nPos, 1 ) ne '*' )
		{
			$rNewQualArray->[$nPos] = $rContig->{base_qualities}[$nQualIndex];
			$nQualIndex++;
		}
		else
		{
			$rNewQualArray->[$nPos] = 0;
		}
	}
	$rContig->{base_qualities} = $rNewQualArray;

	#yes, I know this is nasty, but it's the only quick fix I could think of
	#I'll create a separate data container later
}

sub CreateUnPaddedConsQualityArray
{
	my ($rContig)     = @_;
	my $rNewQualArray = [];
	my $nQualIndex    = 0;
	for ( my $nPos = 0 ; $nPos < length( $rContig->{consensus} ) ; $nPos++ )
	{
		if ( substr( $rContig->{consensus}, $nPos, 1 ) ne '*' )
		{
			$rNewQualArray->[$nQualIndex] = $rContig->{base_qualities}[$nPos];
			$nQualIndex++;
		}
	}
	$rContig->{base_qualities} = $rNewQualArray;

	#yes, I know this is nasty, but it's the only quick fix I could think of
	#I'll create a separate data container later

}

sub SetConsensusToProperCapitals
{
	my ($rContig) = @_;
	for ( my $nPos = 0 ; $nPos < length( $rContig->{consensus} ) ; $nPos++ )
	{
		if ( $rContig->{base_qualities}[$nPos] <= 25 )
		{
			substr( $rContig->{consensus}, $nPos, 1 ) = lc( substr( $rContig->{consensus}, $nPos, 1 ) );
		}
		else
		{
			substr( $rContig->{consensus}, $nPos, 1 ) = uc( substr( $rContig->{consensus}, $nPos, 1 ) );
		}
	}
}

sub SetSequenceToProperCapitals
{
	my ($rRead) = @_;
	if ( !defined( $rRead->{phd} ) || !defined( $rRead->{phd}{sequence} ) || !defined( $rRead->{phd}{sequence}{quality} ) ) { return 1; }
	for ( my $nPos = 0 ; $nPos < length( $rRead->{sequence} ) ; $nPos++ )
	{
		if ( $rRead->{phd}{sequence}{quality}[$nPos] <= 25 )
		{
			substr( $rRead->{sequence}, $nPos, 1 ) = lc( substr( $rRead->{sequence}, $nPos, 1 ) );
		}
		else
		{
			substr( $rRead->{sequence}, $nPos, 1 ) = uc( substr( $rRead->{sequence}, $nPos, 1 ) )
			  unless ( substr( $rRead->{sequence}, $nPos, 1 ) =~ /x|n/ );
		}
	}
}

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

sub numerically { $a <=> $b; }

#returns the read substring index from the consensus position
sub nGetReadIndexFromConsPos
{
	my ( $rRead, $nConsPos ) = @_;

	#the read index is nConsPos + offset, in this case, the offset is
	#(consstart - readstart)
	return ( $nConsPos + ( 1 - $rRead->{read_position}{position} ) );    #removed -1 at the end
}

sub nGetConsPosFromReadIndex
{
	my ( $rRead, $nReadPos ) = @_;
	return ( ( $rRead->{read_position}{position} - 1 ) + $nReadPos );    #removed +1 at the end
}

sub nGetAlignStart
{
	my ($rRead) = @_;
	return $rRead->{read_position}{position};
}

sub nGetAlignEnd
{
	my ($rRead) = @_;

	#ok, this is a bit confusing, read position is in consensus units, and we want the
	#last element in consensus units, so we add base_count-1, not base_count
	return $rRead->{read_position}{position} + $rRead->{padded_base_count} - 1;
}

sub _bIsDyePrimerNotDyeTerminator
{
	my ($rRead) = @_;

	if ( defined( $rRead->{description}{CHEM} ) )
	{
		if ( $rRead->{description}{CHEM} =~ /prim/ )
		{
			return 1;
		}
		elsif ( $rRead->{description}{CHEM} =~ /term/ )
		{
			return 0;
		}
	}

	# if reached here, there are 2 possibilities:
	# Either the CHEM field was not set
	# or, if it was, it was set to something other than
	# prim or term (such as unknown or something unrecognizeable )

	# Thus use the filename extension to detect it.
	my $sFileName = $rRead->{description}{CHROMAT_FILE};

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

sub _bComp
{
	my ($rRead) = @_;

	if ( $rRead->{read_position}{u_or_c} =~ /u|U/ )
	{
		return 0;
	}
	else
	{
		return 1;
	}
}

sub _nGetStrandAndChemistry
{
	my ($rRead) = @_;
	my $bPrimerNotTerminator = _bIsDyePrimerNotDyeTerminator($rRead);

	my $bBottomStrand = _bComp($rRead);

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

sub nFindIndexOfMatchOrPredecessor
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
	if ( defined( $rRead->{phd}{template} ) )    #check phd template first
	{

		if ( $rRead->{phd}{template} =~ /bac|cos/ )
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
sub bAreTheseReadsFromTheSameTemplate
{

	my ( $rRead1, $rRead2 ) = @_;

	if (   _bIsWholeCloneTemplateNotSubcloneTemplate($rRead1)
		|| _bIsWholeCloneTemplateNotSubcloneTemplate($rRead2) )
	{
		return 0;
	}

	# if reached here, neither are whole clone templates

	if ( $rRead1->{phd}{template} eq $rRead2->{phd}{template} )
	{
		return 1;
	}
	else
	{
		return 0;
	}
}

sub _bIsInHighQualitySegmentOfRead
{
	my ( $rRead, $nConsPos ) = @_;
	if ( $rRead->{qual_clip_start} == -1 && $rRead->{qual_clip_end} == -1 )
	{
		return 0;
	}
	else
	{
		if ( $rRead->{qual_clip_start} <= nGetReadIndexFromConsPos( $rRead, $nConsPos )
			&& nGetReadIndexFromConsPos( $rRead, $nConsPos ) <= $rRead->{qual_clip_end} )
		{
			return 1;
		}
		else
		{
			return 0;
		}
	}
}

sub recalculateConsensusQualitiesNoChange
{
	my ( $nRegionConsPosLeft, $nRegionConsPosRight, $rNewQualitiesArray, $rReadArray, $rContig ) = @_;

	$#{$rNewQualitiesArray} = -1;
	my @aReadBoundaries;

	my $nGuessOfNumberOfBoundariesNeeded = ( ( $nRegionConsPosRight - $nRegionConsPosLeft ) * ( $#$rReadArray + 1 ) ) / $rContig->{base_count};

	# the division site will be .5 base to the left of the number
	# put into this array

	push( @aReadBoundaries, $nRegionConsPosLeft );
	push( @aReadBoundaries, ( $nRegionConsPosRight + 1 ) );

	my $nRead;
	for ( $nRead = 0 ; $nRead < @$rReadArray ; $nRead++ )
	{
		my $rRead = $$rReadArray[$nRead];

		#if ( pLocFrag->bIsWholeReadUnaligned() ) continue;
		#if ( pLocFrag->bIsAFakeRead() ) continue;

		my $nIntersectLeft;
		my $nIntersectRight;

		if (
			!bIntersect(
				$nRegionConsPosLeft, $nRegionConsPosRight,
				nGetConsPosFromReadIndex( $rRead, $rRead->{align_clip_start} ),
				nGetConsPosFromReadIndex( $rRead, $rRead->{align_clip_end} ),
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
	@aReadBoundaries = sort numerically @aReadBoundaries;

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

	#$#aConsensusSegments =  $#@aReadBoundaries;
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

		#      if ( pLocFrag->bIsWholeReadUnaligned() ) continue;
		#      if ( pLocFrag->bIsAFakeRead() ) continue;

		my $nIntersectLeft;
		my $nIntersectRight;

		if (
			!bIntersect(
				$nRegionConsPosLeft, $nRegionConsPosRight,
				nGetConsPosFromReadIndex( $rRead, $rRead->{align_clip_start} ),
				nGetConsPosFromReadIndex( $rRead, $rRead->{align_clip_end} ),
				\$nIntersectLeft, \$nIntersectRight
			)
		  )
		{
			next;
		}

		my $nFirstConsSeg = nFindIndexOfMatchOrPredecessor( \@aReadBoundaries, $nIntersectLeft );

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

			if ( bIntervalsIntersect( $nIntersectLeft, $nIntersectRight, $rConsSeg->{cons_pos_left}, $rConsSeg->{cons_pos_right} ) )
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
				if ( !( substr( $rContig->{consensus}, $nConsPos3 - 1, 1 ) eq '*' ) )    #if base at position $nConsPos3 is a pad then...
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
			for ( $nConsPos3 = $nConsPos + 1 ; $nConsPos3 < $rContig->{base_count} ; $nConsPos3++ )
			{
				if ( !( substr( $rContig->{consensus}, $nConsPos3 - 1, 1 ) eq '*' ) )    #if base at position $nConsPos3 is a pad then...
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
						nGetConsPosFromReadIndex( $rRead, $rRead->{align_clip_start} ) <= $nWindowLeft
						&& $nWindowRight <= nGetConsPosFromReadIndex( $rRead, $rRead->{align_clip_end} )
					)
				  )
				{
					next;
				}

				my $bAgreeSoFar = 1;

				for ( my $nConsPos2 = $nWindowLeft ; $nConsPos2 <= $nWindowRight && $bAgreeSoFar ; $nConsPos2++ )
				{
					if (
						lc( substr( $rRead->{sequence}, ( nGetReadIndexFromConsPos( $rRead, $nConsPos2 - 1 ) ), 1 ) ) ne
						lc( substr( $rContig->{consensus}, ( $nConsPos2 - 1 ), 1 ) ) )
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

					my $nQual = $rRead->{phd}{sequence}{quality}[ nGetReadIndexFromConsPos( $rRead, $nConsPos - 1 ) ];

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
							if ( $aTemplates[$nIndex] eq $rRead->{phd}{template} )
							{
								$bMatchFound = 1;
							}
						}
						push( @aTemplates, $rRead->{phd}{template} )
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
							  if ( bAreTheseReadsFromTheSameTemplate( $rReadStrandChem1, $rReadStrandChem2 ) );

							# if reach here, the reads are from different
							# templates

							my $nQual =
							  $rReadStrandChem1->{phd}{sequence}{quality}[ nGetReadIndexFromConsPos( $rReadStrandChem1, $nConsPos - 1 ) ] +
							  $rReadStrandChem2->{phd}{sequence}{quality}[ nGetReadIndexFromConsPos( $rReadStrandChem2, $nConsPos - 1 ) ];

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

					if ( _bIsInHighQualitySegmentOfRead( $rRead, $nConsPos ) )
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
								if ( $aAlignedTemplates[$nIndex] eq $rRead->{phd}{template} )
								{
									$bMatchFound = 1;
								}
							}
						}
						if ( !$bMatchFound )
						{
							push( @aAlignedTemplates, $rRead->{phd}{template} );
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
			if ( @{$rNewQualitiesArray} == 398 )
			{
				my $temp = 1;
			}
			push( @{$rNewQualitiesArray}, $nQualBest );

		}    # nConsPos
	}    # for( int nSegment = 0; nSegment < aConsensusSegments.length();

}

sub recalculateConsensusQualitiesAndChange
{
	my ( $rContig, $nConsPosLeft, $nConsPosRight, $rReadArray ) = @_;

	my @aNewQualities;

	recalculateConsensusQualitiesNoChange( $nConsPosLeft, $nConsPosRight, \@aNewQualities, $rReadArray, $rContig );

	for ( my $nConsPos = $nConsPosLeft ; $nConsPos <= $nConsPosRight ; $nConsPos++ )
	{
		my $qualNew = $aNewQualities[ $nConsPos - $nConsPosLeft ];

		assert( QUALITY_MIN <= $qualNew );
		assert( $qualNew <= QUALITY_MAX );
		$rContig->{base_qualities}[ $nConsPos - 1 ] = $qualNew;
	}
}

sub bGetBSArrayStructureOk
{
	my ( $bCheckThatBaseSegmentsGoFromBeginningToEndOfConsensus, $rBSArray, $rContig ) = @_;

	for ( my $nSeg = 0 ; $nSeg < @$rBSArray ; $nSeg++ )
	{
		if ( $rBSArray->[$nSeg]{start_pos} > $rBSArray->[$nSeg]{end_pos} )
		{
			print "base segment", $nSeg, "has start_pos = ", $rBSArray->[$nSeg]{start_pos}, "and end_pos = ", $rBSArray->[$nSeg]{end_pos},
			  " but start should be <= end\n";

			#need to fix later so that error message prints out base segments corresponding
			#contig
			print "base segment for read ", $rBSArray->[$nSeg]{read_name}, " in contig ";
			return 0;
		}

		# check that the base segment for a read lies within the read

		#first, get the BS's read;
		my $rRead;
		$rRead = $rBSArray->[$nSeg]{read};

		my $nConsPosStart = $rBSArray->[$nSeg]{start_pos};
		my $nConsPosEnd   = $rBSArray->[$nSeg]{end_pos};

		if ( !( nGetAlignStart($rRead) <= $nConsPosStart && $nConsPosStart <= $nConsPosEnd && $nConsPosEnd <= nGetAlignEnd($rRead) ) )
		{
			print "base segment from padded cons pos ", $nConsPosStart, " to ", $nConsPosEnd, " is not within read ", $rRead->{name},
			  " which lies within padded cons pos ", nGetAlignStart($rRead), " to ", nGetAlignEnd($rRead), "\n";
			print "$rRead->{sequence} $rRead->{padded_base_count}\n";
			print "$rContig->{name}\n";
			return 0;
		}

		# skip test the first time
		if ( $nSeg != 0 )
		{
			if ( $rBSArray->[$nSeg]{start_pos} <= $rBSArray->[ $nSeg - 1 ]{end_pos} )
			{
				print "not strictly increasing at ", $nSeg, "\n";
				for ( my $j = ( $nSeg - 1 ) ; ( $j < $nSeg + 2 ) && ( $j < @$rBSArray ) ; $j++ )
				{
					print "index ", $j;
				}

				my $nConsPosStart = $rBSArray->[ $nSeg - 1 ]{start_pos};
				my $nConsPosEnd   = $rBSArray->[ $nSeg + 1 ]{end_post};

				print "base from ", $nConsPosStart, " to ", $nConsPosEnd, "\n";

				print substr( $rContig->{consensus}, $nConsPosStart - 1, ( $nConsPosEnd - $nConsPosStart + 1 ) ), "\n";

				return 0;
			}
		}
	}    # loop through all base segments

	# check the first and last base segments to be sure that the
	# entire consensus is covered with base segments

	if ($bCheckThatBaseSegmentsGoFromBeginningToEndOfConsensus)
	{
		if ( $rBSArray->[0]->{start_pos} != 1 )
		{
			print "Base segment 0 of contig ", $rContig->{name}, " is not at position 1--it should be\n";
			return 0;
		}

		#we need to make sure that BSArray and the contig end at the same pos
		#we use base_count for the contig's end position.  Since the contig
		#index starts at 1, base_count and the end position are the same number
		if ( $rBSArray->[ ( $#{$rBSArray} ) ]->{end_pos} != $rContig->{base_count} )
		{
			print "In Contig ", $rContig->{name}, " last base segment (", ( $#{$rBSArray} ), ") should end on the last padded consensus base (",
			  $rContig->{end_pos}, ") but instead ends on ", $rBSArray->[ $#{$rBSArray} ]->{end_pods}, "\n";

			return 0;
		}
	}

	#if you got here, the array is ok
	return 1;
}

sub max
{
	my ( $x, $y ) = @_;
	if ( $x >= $y )
	{
		return $x;
	}
	return $y;
}

sub min
{
	my ( $x, $y ) = @_;
	if ( $x <= $y )
	{
		return $x;
	}
	return $y;
}

sub nGetBaseSegmentIndexByPos
{
	my ( $rBSarray, $nSeqPos ) = @_;
	my $nIndex;
	my $bFound = 0;
	for ( $nIndex = 0 ; $nIndex < @$rBSarray ; $nIndex++ )
	{
		if (   ( $rBSarray->[$nIndex]{start_pos} <= $nSeqPos )
			&& ( $nSeqPos <= $rBSarray->[$nIndex]{end_pos} ) )
		{
			$bFound = 1;
			last;
		}
	}
	if ($bFound)
	{
		return $nIndex;
	}
	else
	{
		return -1;
	}
}

sub bCheckContiguous
{

	my ( $rReads, $nLeftEndOfTear, $nRightEndOfTear ) = @_;

	my $nPos;
	my @aCoveredByReadInTornRegion;
	for ( $nPos = $nLeftEndOfTear ; $nPos <= $nRightEndOfTear ; $nPos++ )
	{
		$aCoveredByReadInTornRegion[ $nPos++ ] = 0;
	}

	for ( my $nRead = 0 ; $nRead < @$rReads ; $nRead++ )
	{
		my $rRead = $rReads->[$nRead];

		# I tried to decide whether to allow a tear when there might
		# only be unaligned sequence holding one of the resulting
		# contigs together.  I decided 'yes', because there would be
		# unaligned at the ends of each read, hence at the end of each
		# resulting contig.  That is normal, and should be allowed.

		my $nReadInTornRegionLeft;
		my $nReadInTornRegionRight;

		assert(
			bIntersect( $nLeftEndOfTear, $nRightEndOfTear, nGetAlignStart($rRead), nGetAlignEnd($rRead), \$nReadInTornRegionLeft, \$nReadInTornRegionRight ) );

		for ( $nPos = $nReadInTornRegionLeft ; $nPos <= $nReadInTornRegionRight ; $nPos++ )
		{
			$aCoveredByReadInTornRegion[$nPos] = 1;
		}
	}

	for ( $nPos = $nLeftEndOfTear ; $nPos <= $nRightEndOfTear ; $nPos++ )
	{
		if ( !$aCoveredByReadInTornRegion[$nPos] )
		{

			# hole in the contig

			printf( "Contig has position %d that would not be covered by any read.\n", $nPos );
			if ( wantarray() )
			{
				return ( 0, $nPos );
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
		return ( 1, $nPos );
	}
	else
	{
		return 1;
	}
}

sub bIntersect
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

sub bIntervalsIntersect
{
	my ( $nL1, $nR1, $nL2, $nR2 ) = @_;
	return ( ( $nL1 <= $nR2 ) && ( $nL2 <= $nR1 ) );
}

1;
