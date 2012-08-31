package Genome::Model::Tools::Pcap::Contig;
our $VERSION = 0.01;

use strict;
use warnings;
use Carp;

use Genome::Model::Tools::Pcap::Tag;
use List::Util qw(min max);
use Storable;
use Genome::Model::Tools::Pcap::Utility;
use base (qw(Genome::Model::Tools::Pcap::SequenceItem));

=pod

=head1 NAME

Contig - Contig Object.

=head1 SYNOPSIS

my $contig = Genome::Model::Tools::Pcap::Contig->new(reads => \%reads, base_segments => \@base_segments, ace_contig => $ace_contig,  contig_tags => \@contig_tags);

=head1 DESCRIPTION

Genome::Model::Tools::Pcap::Contig mainly acts a container for the reads and sequence data that are normally associated with contig.
It inherits some useful functionality from Genome::Model::Tools::Pcap::SequenceItem.

=head1 METHODS

=cut

=cut
Contig Data Structure
Contig:
	sequence (bases and quality)
	children 
	get_child
	padded_base_count
	base_count
	read_count
	complemented
	tags
=cut

sub new 
{
    croak("__PACKAGE__:new:no class given, quitting") if @_ < 1;
	my ($caller, %params) = @_; 
    my $caller_is_obj = ref($caller);
    my $class = $caller_is_obj || $caller;
    my $self = $class->SUPER::new(%params);
    
	my $ace_contig;    
    
    #eddie suggests to use delete
	$self->{loaded} = 1;
	if(exists $params{callbacks})
	{
		$self->callbacks($params{callbacks});
		$self->{loaded} = 0;
	}
	
	$self->process_params(%params);
	
        
    return $self;
}

sub process_params
{
	my ($self, %params) = @_;
	if(exists $params{reads})
	{
		$self->children ($params{reads});
	}
	if(exists $params{base_segments})
	{
		$self->base_segments ($params{base_segments});
	}
	if(exists $params{contig_tags})
	{
		$self->tags ( $params{contig_tags});
	}
    if(exists $params{ace_contig})
    {
        $self->ace_contig ($params{ace_contig});
    }
}

=pod

=head1 ace_contig

$contig->ace_contig($ace_contig_hash);
    
This is a getter/setter that takes a contig hash that is created by Genome::Model::Tools::Pcap::Ace::Reader.  It will also return an ace contig hash using the data that is contained in the contig.

=cut

sub ace_contig
{
    my ($self, $ace_contig) = @_;

    if(@_ > 1)
    {
        $self->padded_base_string( $ace_contig->{consensus} );
		$self->complemented($ace_contig->{u_or_c} =~ /c/i or 0);
        $self->unpadded_base_quality ( $ace_contig->{base_qualities}); 
        $self->name ($ace_contig->{name});
        $self->{type} = "contig";   
    }
    return { type => "contig",
             name => $self->{name},
             base_count => $self->length,
             read_count => scalar (keys %{$self->reads}),
             base_seg_count => scalar (@{$self->base_segments}),
             u_or_c => ( $self->complemented ? "C" : "U" ),
             consensus => $self->padded_base_string,
             base_qualities => $self->unpadded_base_quality
           };             
}

=pod

=head1 reads

$contig->reads(\%reads);
    
This is a getter/setter that will either get or set the Contig's hash of reads.

=cut

sub reads
{
	my ($self, $reads) = @_;

	$self->children( $reads ) if (@_ > 1);	
	 
	return $self->children;
}

=pod

=head1 base_segments

$contig->base_segments(\@base_segments);
    
This is a getter/setter that will either get or set the Contig's array of base_segments.

=cut

sub base_segments
{
	my ($self, $value) = @_;
	my ($name) = (caller(0))[3] =~ /.+::(.+)/;
	
    if(@_>1)
    {   
       return $self->check_and_load_data($name, $value);
    }
    return $self->check_and_load_data($name);
	
}

=pod

=head1 read_count

$contig->read_count;
    
This returns the count of reads.

=cut

sub read_count
{
	my ($self, $value) = @_;
	
	if(defined $self->{children})
    {
        return scalar keys %{$self->{children}};
    }    
    my ($name) = (caller(0))[3] =~ /.+::(.+)/;
    if(@_>1)
    {   
       return $self->check_and_load_data($name, $value);
    }
    return $self->check_and_load_data($name);
}

=pod

=head1 base_count

$contig->base_count;
    
This returns the number of unpadded bases.

=cut

sub base_count
{
	my ($self, $value) = @_;
    if($self->already_loaded("padded_base_string"))
    {
	   return length $self->padded_base_string;
    }
    my ($name) = (caller(0))[3] =~ /.+::(.+)/;
    if(@_>1)
    {   
       return $self->check_and_load_data($name, $value);
    }
    return $self->check_and_load_data($name);
}

=pod

=head1 base_segment_count

my $bs_count = $contig->base_segment_count;
    
This returns the number of base segments.

=cut

sub base_segment_count
{
	my ($self, $value) = @_;
	if(defined $self->{base_segments})
    {
	   return scalar @{$self->base_segments};
    }
    my ($name) = (caller(0))[3] =~ /.+::(.+)/;
    if(@_>1)
    {   
       return $self->check_and_load_data($name, $value);
    }
    return $self->check_and_load_data($name);
    
}

=pod

=head1 complemented

my $iscomplemented = $contig->complemented($iscomplemented);

=cut

sub complemented
{
	my ($self, $value) = @_;
	
    my ($name) = (caller(0))[3] =~ /.+::(.+)/;
    if(@_>1)
    {   
       return $self->check_and_load_data($name, $value);
    }
    return $self->check_and_load_data($name);
}

sub is_bs_array_structure_ok
{
	my ( $self, $bCheckThatBaseSegmentsGoFromBeginningToEndOfConsensus, $rBSArray ) = @_;
	if(!defined $rBSArray) { $rBSArray = $self->base_segments;}
	my $reads = $self->reads;
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
		$rRead = $reads->{$rBSArray->[$nSeg]{read_name}};
		my $nConsPosStart = $rBSArray->[$nSeg]{start_pos};
		my $nConsPosEnd   = $rBSArray->[$nSeg]{end_pos};

		if ( !( $rRead->align_start <= $nConsPosStart && $nConsPosStart <= $nConsPosEnd && $nConsPosEnd <= $rRead->align_end ) )
		{
			print "base segment from padded cons pos ", $nConsPosStart, " to ", $nConsPosEnd, " is not within read ", $rRead->{name},
			  " which lies within padded cons pos ", $rRead->align_start, " to ", $rRead->align_end, "\n";
			print "$rRead->padded_base_string() $rRead->base_count()\n";
			print "$self->name()\n";
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
				my $nConsPosEnd   = $rBSArray->[ $nSeg + 1 ]{end_pos};

				print "base from ", $nConsPosStart, " to ", $nConsPosEnd, "\n";

				print substr( $self->padded_base_string, $nConsPosStart - 1, ( $nConsPosEnd - $nConsPosStart + 1 ) ), "\n";

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
			print "Base segment 0 of contig ", $self->name, " is not at position 1--it should be\n";
			return 0;
		}

		#we need to make sure that BSArray and the contig end at the same pos
		#we use base_count for the contig's end position.  Since the contig
		#index starts at 1, base_count and the end position are the same number
		if ( $rBSArray->[ ( $#{$rBSArray} ) ]->{end_pos} != $self->base_count )
		{
			print "In Contig ", $self->name, " last base segment (", ( $#{$rBSArray} ), ") should end on the last padded consensus base (",
			  $self->base_count, ") but instead ends on ", $rBSArray->[ $#{$rBSArray} ]->{end_pos}, "\n";

			return 0;
		}
	}

	#if you got here, the array is ok
	return 1;
}


sub children
{
    my ($self, $value) = @_;
    
    my ($name) = (caller(0))[3] =~ /.+::(.+)/;
    if(@_>1)
    {   
       return $self->check_and_load_data($name, $value);
    }
    return $self->check_and_load_data($name);
}

=pod

=head1 copy_base_segment

my $base_segment_copy = $contig->copy_base_segment($base_segment);
    
This function will copy a base_segment.

=cut

sub copy_base_segment
{
	my ($self, $base_segment) = @_;
	return dclone($base_segment);	
}

=pod

=head1 copy_tag

my $tag_copy = $contig->copy_tag($contig_tag);
    
This copies a contig tag.

=cut

=pod

=head1 copy

my $contig_copy = $contig->copy($contig);
    
This creates a deep copy of the contig.

=cut

sub calculate_base_segments
{
	my ($self, $start_pos, $end_pos,$best_quality_reads) = @_;
	if(!defined $best_quality_reads)
	{
		($best_quality_reads) = $self->get_best_quality_reads($start_pos,$end_pos);	
	}
	my @bs;
	for ( my $pos = $start_pos ; #we use cons_start_new since we need to account for the new length
		$pos <= $end_pos ;  #of the merge region
		$pos++ )
	{
		my $base_seg = {
			type      => 'base_segment',
			start_pos => $pos ,
			end_pos   => $pos ,
			read_name => $best_quality_reads->[$pos]->name,
		};
		push( @bs, $base_seg );
	}
	return \@bs;
}

sub calculate_consensus
{
	my ($self, $start_pos, $end_pos,$overlap_reads) = @_;
	my $consensus;
	my @quality;
	my ($best_quality_reads, $best_quality) = $self->get_best_quality_reads($start_pos,$end_pos,$overlap_reads);
	
	for ( my $pos = $start_pos; $pos <= $end_pos; $pos++ )
	{
		my $read = $best_quality_reads->[$pos];
		$consensus .= substr( $read->padded_base_string, $read->get_child_position_from_parent_position($pos-1), 1 ); 
		$quality[$pos-1] = $best_quality->[$pos];
		#extend align region for reads, I'm not sure if this belongs here or not
		if ( $pos < $read->get_parent_position_from_child_position($read->align_clip_start))
		{
			$read->align_clip_start($read->get_child_position_from_parent_position( $pos ));
		}
		elsif($pos > $read->get_parent_position_from_child_position( $read->align_clip_end ))
		{
			$read->align_clip_end($read->get_child_position_from_parent_position( $pos ));
		}
	}
	$self->padded_base_string($consensus);
	$self->padded_base_quality(\@quality);
	return ($best_quality_reads, $best_quality);
}


sub recalculate_consensus
{
	my ($self, $start_pos, $end_pos,$overlap_reads) = @_;
	my $consensus = $self->padded_base_string;
	my @quality;
	eval { @quality = @{$self->padded_base_quality}; };
	if(@quality == 0)
	{
		for(my $i=0;$i<$self->length;$i++){$quality[$i]=0;}
	}
	my ($best_quality_reads, $best_quality) = $self->get_best_quality_reads($start_pos,$end_pos,$overlap_reads);
	
	for ( my $pos = $start_pos; $pos <= $end_pos; $pos++ )
	{
		my $read = $best_quality_reads->[$pos];
		substr($consensus, $pos-1, 1) = substr( $read->padded_base_string, $read->get_child_position_from_parent_position($pos-1), 1 ); 
		$quality[$pos-1] = $best_quality->[$pos];
		#extend align region for reads, I'm not sure if this belongs here or not
		if ( $pos < $read->get_parent_position_from_child_position($read->align_clip_start))
		{
			$read->align_clip_start($read->get_child_position_from_parent_position( $pos ));
		}
		elsif($pos > $read->get_parent_position_from_child_position( $read->align_clip_end ))
		{
			$read->align_clip_end($read->get_child_position_from_parent_position( $pos ));
		}
	}
	$self->padded_base_string($consensus);
	$self->padded_base_quality(\@quality);
	return ($best_quality_reads, $best_quality);
}

#this was dragged over from Utility.pm to remove the dependency on this library.
use constant QUALITY_MIN               => 0;
use constant QUALITY_LOW_EDITED        => 98;
use constant QUALITY_HIGH_EDITED       => 99;
sub nNormalQualityFrom9899Quality
{
	my ($nQuality) = @_;
    if($nQuality =~ /\D+/)
    {
        $nQuality = 0;
    }
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

sub get_best_quality_reads
{
	my ($self, $start_pos, $end_pos,$overlap_reads) = @_;
	my @aBestQualityRead;
	my @aBestQuality;
	$#aBestQualityRead = -1;
	$#aBestQuality = -1;
	$#aBestQualityRead = $end_pos- $start_pos+1;
	$#aBestQuality = $end_pos-$start_pos+1;

	if(!defined $overlap_reads)
	{
		$overlap_reads = $self->reads;	
	}
	foreach my $read (values %{$overlap_reads})
	{
		my $left_intersect;
		my $right_intersect;
		
		#assert(
		next if (!
		calculate_intersection(
			$start_pos,
			$end_pos,
			$read->align_start,
			$read->align_end,
			\$left_intersect,
			\$right_intersect
		));
		#);
		my $qual_clip_left= $read->get_parent_position_from_child_position($read->qual_clip_start);
		my $qual_clip_right= $read->get_parent_position_from_child_position($read->qual_clip_end);

		for ( my $nConsPos = $left_intersect ;
			$nConsPos <= $right_intersect ;
			$nConsPos++ )
		{
			if ( !$aBestQualityRead[$nConsPos] )
			{
				$aBestQualityRead[$nConsPos] = $read;

				if ( lc(substr( $read->padded_base_string, $read->get_child_position_from_parent_position( $nConsPos-1), 1 )) eq 'x' )
				{
					$aBestQuality[$nConsPos] = 0;
				}
				else
				{					
					$aBestQuality[$nConsPos] =
					
					nNormalQualityFrom9899Quality(
						$read->padded_base_quality->
						[ $read->get_child_position_from_parent_position( $nConsPos-1 ) ] );
				}
			}
			else
			{
				my $q;
				if ( lc(substr( $read->padded_base_string, $read->get_child_position_from_parent_position( $nConsPos-1), 1 )) eq 'x' )
				{
					$q = 0;
				}
				else
				{					
					$q =
					nNormalQualityFrom9899Quality(
						$read->padded_base_quality->
						[ $read->get_child_position_from_parent_position( $nConsPos-1 ) ] );
				}

				if ( $q > $aBestQuality[$nConsPos] )
				{
					$aBestQualityRead[$nConsPos] = $read;
					$aBestQuality[$nConsPos]     = $q;
				}			
			}
		}
	}
	return (\@aBestQualityRead,\@aBestQuality);
}

sub get_overlapping_reads
{
	my ($self, $start_position, $end_position, $reads) = @_;
	my %overlap_reads;
	$reads = values %{$self->reads} if !defined $reads;
	foreach my $read ( values %{$reads})
	{	
		$overlap_reads{$read->name} = $read if (
			 intervals_intersect( $start_position,
								  $end_position,
								  $read->align_start,
								  $read->align_end
			 )
		  );		
	}
	return \%overlap_reads;
}

sub loaded
{
	my ($self, $loaded) = @_;
	$self->{loaded} = $loaded if (@_>1);
	return $self->{loaded};
}

sub _get_extents_for_reads
{
	my ($self, $reads) = @_;
	my $init = 1;
	my ($cons_start, $cons_end);
	foreach my $read (@{$reads})
	{
		my $read_start = $read->get_parent_position_from_child_position($read->align_clip_start);
		my $read_end = $read->get_parent_position_from_child_position($read->align_clip_end);
		if($init)
		{
			$init = 0;
			($cons_start, $cons_end) = ($read_start,$read_end); 
		}
		$cons_start = min($read_start,$cons_start);
		$cons_end = max($read_end,$cons_end);
	}
	return ($cons_start, $cons_end);



}
sub get_contiguous_reads
{
	my ($self, $reads, $tags) = @_;
	my @reads = sort {$a->position <=> $b->position} values %{$reads};
	my @coverage;
	my ($cons_start,$cons_end) = $self->_get_extents_for_reads(\@reads);
	#get offset
	my $offset = 1-$cons_start;
	#move the reads over by $offset
	foreach (@reads) { $_->position($_->position+$offset);} 
	$cons_start+=$offset;
	$cons_end+=$offset;
	foreach my $read (@reads)
	{
		my $read_start = $read->get_parent_position_from_child_position($read->align_clip_start);
		my $read_end = $read->get_parent_position_from_child_position($read->align_clip_end);
		for(my $i = $read_start;$i<=$read_end;$i++)
		{
			if (!defined $coverage[$i])
			{
				$coverage[$i] = [$read];
			}
			else
			{
				push @{$coverage[$i]}, $read;
			}
		}
	}
	my @cont_regions;
	my $region_index = 0;
	for(my $i = $cons_start;$i<=$cons_end;$i++)
	{
		if(defined $coverage[$i])
		{
			if(!defined $cont_regions[$region_index])
			{
				$cont_regions[$region_index]->{start_pos} = $i;
				$cont_regions[$region_index]->{length} = 0;				
			}
			$cont_regions[$region_index]->{end_pos} = $i;
			$cont_regions[$region_index]->{length}++;
			
			foreach my $read (@{$coverage[$i]})
			{
				$cont_regions[$region_index]->{reads}{$read->name} = $read;				
			}		
		}
		elsif(defined $coverage[$i-1])
		{
			$region_index++;
		}			
	}
#	return @cont_regions;
	my $longest_region = $cont_regions[0];
	if(@cont_regions>1)
	{
		foreach my $region (@cont_regions)
		{
			$longest_region = $region if($region->{length}>$longest_region->{length});			
		}
		#go ahead and figure out offset
				
		my $offset = 1-$longest_region->{start_pos};
		
		#move the reads over by $offset
		foreach my $read (values %{$longest_region->{reads}})
		{ $read->position($read->position+$offset);}
		$longest_region->{offset} = $offset;
		#$longest_region->{start_pos};
		#$longest_region->{end_pos};
		if((scalar @{$tags})>0)
		{
			$longest_region->{tags} = $self->clip_and_filter_tags($longest_region,$tags);
		}
	}
	return $longest_region;
}

sub clip_and_filter_tags
{
	my ($self,$region, $tags) = @_;
	my ($new_start,$new_end);
	my @new_tags;
	foreach my $tag (@$tags)
	{
		if(calculate_intersection($tag->start,$tag->stop,$region->{start_pos},$region->{end_pos},\$new_start,\$new_end))
		{
			$tag->start($new_start+$region->{offset});
			$tag->stop($new_end+$region->{offset});
			push @new_tags,$tag;
		}	
	}
	return \@new_tags;
}


sub get_bs_index_at_position
{
	my ( $self, $nSeqPos, $bs ) = @_;
	
	$bs = $self->base_segments if( ! defined ($bs));
	my $nIndex;
	my $bFound = 0;
	for ( $nIndex = 0 ; $nIndex < @$bs ; $nIndex++ )
	{
		if (   ( $bs->[$nIndex]{start_pos} <= $nSeqPos )
			&& ( $nSeqPos <= $bs->[$nIndex]{end_pos} ) )
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

sub extend
{
	my ($self, $side, $po) = @_;
	if($side eq "left")
	{
		my %reads = %{$self->reads};
		my @reads;
		foreach my $read (values %reads)
		{
			push @reads, $read if($read->align_start < 1 && $read->qual_clip_start != -1 && $read->align_clip_start != -1);
		}
		
		foreach my $read (@reads)		
		{
			my $phd =$po->get_phd($read->phd_file);
			my @quality = @{$phd->unpadded_base_quality};
			@quality = reverse @quality if ($read->complemented);
			$read->unpadded_base_quality( \@quality);
			my $qual_array = $read->padded_base_quality;
			my $start_cons_in_read = $read->get_child_position_from_parent_position(1);
			my $sum = 0;
			for(my $i = 1;$i<$start_cons_in_read;$i++)
			{
				next if($qual_array->[$i] eq '*');
				$sum += $qual_array->[$i];
			}
			$read->{sum} = $sum;
		}
		@reads = sort { $a->{sum} <=> $b->{sum} } @reads;
		my $extend_read = pop @reads;
		foreach my $read (@reads)
		{
			$read->padded_base_quality([]);
		}		
		
		$extend_read->qual_clip_start(1);
		$extend_read->align_clip_start(1);
		
		my $length = 1 - $extend_read->align_start;
		my $newbs = { type => 'base_segment', 
					  start_pos => 1,
					  end_pos => $length,
					  read_name => $extend_read->name };
		my @bs = @{$self->base_segments};
		foreach my $bs (@bs)
		{
			$bs->{start_pos} += $length;
			$bs->{end_pos} += $length;		
		}
		unshift @bs, $newbs;
		$self->base_segments(\@bs);
		my $consensus = substr($extend_read->padded_base_string,0,$length).$self->padded_base_string;
		my @temp_quality = @{$extend_read->padded_base_quality};
		@temp_quality = splice(@temp_quality,0,$length);
		push @temp_quality,@{$self->padded_base_quality};
		$self->padded_base_string($consensus);
		$self->padded_base_quality(\@temp_quality);
		
		foreach my $read (values %reads)
		{
			$read->position($read->position + $length);
		}
		
		my @tags = @{$self->tags};
		foreach my $tag (@tags)
		{
			$tag->start($tag->start+$length);
			$tag->stop($tag->stop+$length);
		}
		$self->tags(\@tags);		
	}
	else
	{
		my %reads = %{$self->reads};
		my @reads;
		foreach my $read (values %reads)
		{
			push @reads, $read if($read->align_end > $self->length && $read->qual_clip_start != -1 && $read->align_clip_start != -1);		
		}
		foreach my $read (@reads)		
		{
			my $phd =$po->get_phd($read->phd_file);
			my @quality = @{$phd->unpadded_base_quality};
			@quality = reverse @quality if ($read->complemented);
			$read->unpadded_base_quality( \@quality);
			my $qual_array = $read->padded_base_quality;
			my $end_cons_in_read = $read->get_child_position_from_parent_position($self->length);
			my $sum = 0;
			for(my $i = $end_cons_in_read;$i<$read->length;$i++)
			{
				next if($qual_array->[$i] eq '*');
				$sum += $qual_array->[$i];
			}
			$read->{sum} = $sum;
		}
		@reads = sort { $a->{sum} <=> $b->{sum} } @reads;
		my $extend_read = pop @reads;
		foreach my $read (@reads)
		{
			$read->padded_base_quality([]);
		}		
		
		my $newbs = { type => 'base_segment', 
					  start_pos => ($self->length + 1),
					  end_pos => $extend_read->align_end,
					  read_name => $extend_read->name };
		my @bs = @{$self->base_segments};
		push @bs, $newbs;
		$self->base_segments(\@bs);
		$extend_read->align_clip_end($extend_read->length);
		$extend_read->qual_clip_end($extend_read->length);
		my $length = $extend_read->align_end - $self->length;
		my $read_offset = $extend_read->length - $length;
		my $consensus = $self->padded_base_string.substr($extend_read->padded_base_string,$read_offset,$length);
		my @temp_quality = @{$extend_read->padded_base_quality};
		@temp_quality = splice(@temp_quality,$read_offset,$length);
		unshift @temp_quality,@{$self->padded_base_quality};
		$self->padded_base_string($consensus);
		$self->padded_base_quality(\@temp_quality);			
	}

}


1;

