package Finishing::Assembly::Contig;

use strict;
use warnings;
no warnings 'reserved';

use base 'Finishing::Assembly::SequencedItem';
use List::Util qw(min max);

use Data::Dumper;
use Finishing::Assembly::Utility;

require Storable;

my %align :name(align:o) :ds(hashref) :empty_ok(1) default({}); #added for compatibility with ContigTools

#- ASSEMBLED READS -#
sub reads
{
    my $self = shift;

    return $self->assembled_reads(@_);
}

sub read_count
{
    return shift->assembled_reads->count;
} 

#- TAGS -#
sub tags
{
    my ($self, $tags) = @_;
	my ($tag_conds, $tag_params);# = @_;
    
	my $method = $self->proxy->get_method('tags',$tags);
    if ( $tag_conds )
    {
        $self->set_unpadded_start_and_stop_in_tags if exists $tag_conds->{unpad_start} 
            or exists $tag_conds->{unpad_stop};
 
        return $method->()->search($tag_conds, $tag_params);
    }
    else
    {
        return $method->();
    }
}



sub set_unpadded_start_and_stop_in_tags
{
    my $self = shift;

    my $ti = $self->get_tag_iterator;
    while ( my $tag = $ti->next )
    {
        $tag->unpad_start( $self->unpad_position_to_pad_position( $tag->start ) );
        $tag->unpad_stop( $self->unpad_position_to_pad_position( $tag->stop ) );
    }
    
    return 1;
}

sub add_tags
{
    my ($self, $new_tags) = @_;

    my $tags = $self->tags;

    push @$tags, $new_tags;

    return $self->tags($tags);
}

sub remove_tags
{
    my $self = shift;

    my $ti = $self->get_tag_iterator(@_);
    
    return $self->tags([ $ti->all ]);
}

sub remove_all_tags
{
    my $self = shift;

    return $self->tags([]);
}

sub remove_duplicate_tags
{
    my $self = shift;

    my $tags = $self->tags;
    return unless @$tags;
    
    my %unique_tags;
    foreach my $tag ( @$tags )
    {
        my $key = join
        (
            '', $tag->parent, $tag->type, $tag->program, 
            $tag->start, $tag->stop, 
            $tag->date, $tag->text, $tag->comment,
        );
        $unique_tags{$key} = $tag;
    }
    
    return $self->tags([ values %unique_tags ]);
}

#- COUNTS -#
sub base_segment_count
{
    my $self = shift;

    return scalar @{ $self->base_segments };
}

#- OPERATIONS -#
sub complement : CUMULATIVE
{	
	my ($self) = @_;

    my $length = $self->length;

	my $reads = $self->assembled_reads;
    
	while ( my $read = $reads->next )
	{
        $read->complement;
        $read->position(  $length - ($read->position + $read->length ) + 2 );
    }
    
    my @bs = reverse @{$self->base_segments};
	foreach my $bs ( @bs )
	{
		my $start = $length - $bs->{stop} + 1;
        my $end = $length - $bs->{start} + 1; 
		$bs->{start} = $start;
		$bs->{stop} = $end;
	}
	$self->base_segments(\@bs);
	
	foreach my $tag ( @{ $self->tags } ) 
	{
        my $start = $length - $tag->stop + 1;
        my $end = $length - $tag->start + 1;
		$tag->start($start);
		$tag->stop($end);	
	}
    
    return 1;
}

sub is_bs_array_structure_ok
{
	my ($self, $check_bs_are_contiguous, $base_segments) = @_;

    $base_segments = $self->base_segments unless $base_segments;
    
    my $bs_count = scalar @$base_segments;
    my $bs_last_pos = $bs_count - 1;
	my %reads = map {$_->name, $_} $self->assembled_reads->all;
	for (my $i = 0; $i < $bs_count; $i++)
	{
        my $bs = $base_segments->[$i];
        my $bs_read_name = $bs->{read_name};
        $self->error_msg("No read name for base segment position $i")
            and return unless $bs_read_name;
        my $bs_start = $bs->{start};
        $self->error_msg("No start for base segment position $i, read name $bs_read_name")
            and return unless $bs_start;
        my $bs_stop = $bs->{stop};
        $self->error_msg("No stop for base segment position $i, read name $bs_read_name")
            and return unless $bs_stop;
        my $read = $reads{$bs_read_name};#$self->get_assembled_read($bs_read_name);
        $self->error_msg
        (
            sprintf
            (
                '%s base segment\'s (#%d, read: %s) start (%d) is greater than its stop (%d)',
                $self->name,
                $i,
                $bs_read_name,
                $bs_start,
                $bs_stop,
            )
        ) and return if $bs_start > $bs_stop;

        $self->error_msg
        (
            sprintf
            (
                "%s base segment\'s (#%d, %d to %d) is not w/in its read (%s, %d to %d, %d)\n%s",
                $self->name,
                $i,
                $bs_start,
                $bs_stop,
                $bs_read_name,
                $read->start,
                $read->stop,
                $read->length,
                $read->base_string,
            )
        ) and return unless $read->align_start <= $bs_start 
            and $bs_start <= $bs_stop 
            and $bs_stop <= $read->align_stop;

		if ( $i ) # skip test the first time
		{
            my $prev_bs = $base_segments->[$i - 1];
            my $prev_start = $prev_bs->{start};
            my $prev_stop = $prev_bs->{stop};
			if ( $bs_start <= $prev_stop )
			{
                my $msg = "Not strictly increasing at $i\n";
				for ( my $j = ( $i - 1 ); ( $j < $i + 2 ) and ( $j < $bs_count ); $j++)
				{
                    $msg .= "index : $j\n";
				}
				my $start = $base_segments->[ $i - 1 ]->{start};
				my $stop   = $base_segments->[ $i + 1 ]->{stop};
                $msg .= "base from $start to $stop\n";
                $msg .= substr
                (
                    $self->base_string, 
                    $start - 1,
                    ( $stop - $start + 1 ) 
                ) . "\n";
                $self->error_msg($msg);
                return;
			}
		}
	}   

	if ( $check_bs_are_contiguous )
	{
		unless ( $base_segments->[0]->{start} == 1 )
		{
            $self->error_msg("Base segment 0 of contig " . $self->name . " is not at position 1");
            return;
		}

		unless ( $base_segments->[$bs_last_pos]->{stop} == $self->length )
		{
            $self->error_msg
            (
                sprintf
                (
                    'In Contig %s the last base segment (%d) should end on the last padded consensus base '.
					"(%d) but instead ends on (%d)\n",
					$self->name,
					$self->length,
					$bs_last_pos,
                    $base_segments->[$bs_last_pos]->{stop},                    
                )
            );
            return;
		}
	}

	return 1;
}

sub calculate_base_segments
{
	my ($self, $start, $stop, $best_quality_reads) = @_;

    $best_quality_reads = $self->get_best_quality_reads($start, $stop) unless $best_quality_reads;

	my @bs;
	for ( my $pos = $start; $pos <= $stop; $pos++ )
	{
        push @bs,  
        {
			start => $pos,
			stop   => $pos,
			read_name => $best_quality_reads->[$pos]->name,
		};
	}

	return \@bs;
}


#################################
# not refactored below here#
sub calculate_consensus
{
	my ($self, $start, $stop, $overlap_reads) = @_;

	my $consensus;
	my @quality;
	my ($best_quality_reads, $best_quality) = $self->get_best_quality_reads($start,$stop,$overlap_reads);
	
	for ( my $pos = $start; $pos <= $stop; $pos++ )
	{
		my $read = $best_quality_reads->[$pos];
		$consensus .= substr( $read->padded_base_string, $read->get_child_position_from_parent_position($pos-1), 1 ); 
		$quality[$pos-1] = $best_quality->[$pos];
		#extend align region for reads, I'm not sure if this belongs here or not
		if ( $pos < $read->get_parent_position_from_child_position($read->align_clip_start))
		{
			$read->align_clip_start($read->get_child_position_from_parent_position( $pos ));
		}
		elsif($pos > $read->get_parent_position_from_child_position( $read->align_clip_stop ))
		{
			$read->align_clip_stop($read->get_child_position_from_parent_position( $pos ));
		}
	}
	$self->padded_base_string($consensus);
	$self->padded_base_quality(\@quality);
	return ($best_quality_reads, $best_quality);
}

sub recalculate_consensus
{
	my ($self, $start, $stop,$overlap_reads) = @_;
	my $consensus = $self->padded_base_string;
	my @quality;
	eval { @quality = @{$self->padded_base_quality}; };
	if(@quality == 0)
	{
		my $length = $self->length;
		for(my $i=0;$i<$length;$i++){$quality[$i]=0;}
	}
	my ($best_quality_reads, $best_quality) = $self->get_best_quality_reads($start,$stop,$overlap_reads);
	#print $self->name."\n";
	for ( my $pos = $start; $pos <= $stop; $pos++ )
	{
	
		my $read = $best_quality_reads->[$pos];
		#print $read->name."\n";
		substr($consensus, $pos-1, 1) = substr( $read->padded_base_string, $read->get_child_position_from_parent_position($pos-1), 1 ); 
		$quality[$pos-1] = $best_quality->[$pos];
		#extend align region for reads, I'm not sure if this belongs here or not
		if ( $pos < $read->get_parent_position_from_child_position($read->align_clip_start))
		{
			$read->align_clip_start($read->get_child_position_from_parent_position( $pos ));
		}
		elsif($pos > $read->get_parent_position_from_child_position( $read->align_clip_stop ))
		{
			$read->align_clip_stop($read->get_child_position_from_parent_position( $pos ));
		}
	}
	$self->padded_base_string($consensus);
	$self->padded_base_quality(\@quality);
	return ($best_quality_reads, $best_quality);
}

sub get_best_quality_reads
{
	my ($self, $start, $stop, $overlap_reads) = @_;

	my @aBestQualityRead;
	my @aBestQuality;
	$#aBestQualityRead = -1;
	$#aBestQuality = -1;
	$#aBestQualityRead = $stop- $start+1;
	$#aBestQuality = $stop-$start+1;

    $overlap_reads = [$self->assembled_reads->all] unless defined $overlap_reads;
	
	foreach my $read ( @{$overlap_reads} )
	{
		my $left_intersect;
		my $right_intersect;
		
		#assert(
		next if (!
		calculate_intersection(
			$start,
			$stop,
			$read->align_start,
			$read->align_stop,
			\$left_intersect,
			\$right_intersect
		));
		#);
		my $qual_clip_left= $read->get_parent_position_from_child_position($read->qual_clip_start);
		my $qual_clip_right= $read->get_parent_position_from_child_position($read->qual_clip_stop);

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
					
					Finishing::Assembly::Utility::nNormalQualityFrom9899Quality(
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
					Finishing::Assembly::Utility::nNormalQualityFrom9899Quality(
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
	my ($self, $start, $stop, $reads) = @_;	

	$reads = [$self->assembled_reads->all] unless $reads;

	my @overlap_reads;
	foreach my $read (@$reads)
	{	
        push @overlap_reads, $read if intervals_intersect
        (
            $start,
            $stop,
            $read->align_start,
            $read->align_stop,
        );		
    }

    return \@overlap_reads;
}

sub _get_extents_for_reads
{
	my ($self, $reads) = @_;

	my $init = 1;
	my ($cons_start, $cons_end);
	foreach my $read (@{$reads})
	{
		my $read_start = $read->get_parent_position_from_child_position($read->align_clip_start);
		my $read_end = $read->get_parent_position_from_child_position($read->align_clip_stop);
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
	my @reads = sort {$a->position <=> $b->position} @{$reads};
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
		my $read_end = $read->get_parent_position_from_child_position($read->align_clip_stop);
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
				push @{$cont_regions[$region_index]->{reads}}, $read;				
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
		foreach my $read (@{$longest_region->{reads}})
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
		if (   ( $bs->[$nIndex]->{start} <= $nSeqPos )
			&& ( $nSeqPos <= $bs->[$nIndex]->{stop} ) )
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
	my ($self, $side, $po, $qual_thresh) = @_;
    
	my @all_reads = $self->assembled_reads->all;
    my %sums = ();
    my %avg  = ();
    $qual_thresh ||= 15;
       
	if($side eq "left")
	{
		my @reads;
		foreach my $read (@all_reads)
		{
			push @reads, $read if($read->align_start < 1 && $read->qual_clip_start != -1 && $read->align_clip_start != -1);
		}
		
		foreach my $read (@reads)		
		{
            #my $phd =$po->get_phd($read->phd_file);
            my $phd = get_phd($read->phd_file, $po);
			my @quality = @{$phd->unpadded_base_quality};
			@quality = reverse @quality if ($read->complemented);
			$read->unpadded_base_quality( \@quality);
			my $qual_array = $read->padded_base_quality;
			my $seq = $read->padded_base_string;
			my @seq_array = split //, $seq;
			@seq_array = reverse @seq_array if($read->complemented);
			my $start_cons_in_read = $read->get_child_position_from_parent_position(1);
			my $sum = 0;
            my $ct  = 0;
			for(my $i = 1;$i<$start_cons_in_read;$i++)
			{
				next if($qual_array->[$i] eq '*');
				if($seq_array[$i] =~ /N|n|X|x/){ $sum = 0;last;}
				$sum += $qual_array->[$i];
                $ct++;
			}
            $sums{$sum} = $read;
            $avg{$sum} = $sum / $ct if $ct;
            #$read->{sum} = $sum;
		}
        my @order = sort {$a <=> $b}keys %sums;
        #@reads = sort { $a->{sum} <=> $b->{sum} } @reads;
        $self->warn_msg("Average base quality of extend read is lower than threshold $qual_thresh, No extend")
            and return unless $avg{$order[$#order]} > $qual_thresh;

		my $extend_read = $sums{pop @order};
		
		$extend_read->qual_clip_start(1);
		$extend_read->align_clip_start(1);
		
		my $length = 1 - $extend_read->align_start;
		my $newbs = { type => 'base_segment', 
					  start => 1,
					  stop => $length,
					  read_name => $extend_read->name };
		my @bs = @{$self->base_segments};
		foreach my $bs (@bs)
		{
			$bs->{start} += $length;
			$bs->{stop} += $length;		
		}
		unshift @bs, $newbs;
		$self->base_segments(\@bs);
		my $consensus = substr($extend_read->padded_base_string,0,$length).$self->padded_base_string;
		my @temp_quality = @{$extend_read->padded_base_quality};
		@temp_quality = splice(@temp_quality,0,$length);
		push @temp_quality,@{$self->padded_base_quality};
		$self->padded_base_string($consensus);
		$self->padded_base_quality(\@temp_quality);
		
		foreach my $read (@all_reads)
		{
			$read->position($read->position + $length);
		}
		
		my @tags = @{$self->tags};

        if (@tags) {
		    foreach my $tag (@tags)
		    {
			    $tag->start($tag->start+$length);
			    $tag->stop($tag->stop+$length);
		    }
		    $self->tags(\@tags);	
        }
	}
	else
	{
		my @reads;
		foreach my $read (@all_reads)
		{
			push @reads, $read if($read->align_stop > $self->length && $read->qual_clip_start != -1 && $read->align_clip_start != -1);		
		}
		foreach my $read (@reads)		
		{
            #my $phd =$po->get_phd($read->phd_file);
            my $phd = get_phd($read->phd_file, $po);
			my @quality = @{$phd->unpadded_base_quality};
			@quality = reverse @quality if ($read->complemented);
			$read->unpadded_base_quality( \@quality);
			my $qual_array = $read->padded_base_quality;
			my $seq = $read->padded_base_string;
			my @seq_array = split //, $seq;
			@seq_array = reverse @seq_array if($read->complemented);
			my $end_cons_in_read = $read->get_child_position_from_parent_position($self->length);
			my $sum = 0;
            my $ct  = 0;
			for(my $i = $end_cons_in_read;$i<$read->length;$i++)
			{
				next if($qual_array->[$i] eq '*');
				if($seq_array[$i] =~ /N|n|X|x/){ $sum = 0;last;}
				$sum += $qual_array->[$i];
                $ct++;
			}
            $sums{$sum} = $read;
            #$read->{sum} = $sum;
            $avg{$sum} = $sum / $ct if $ct;
		}
        my @order = sort {$a <=> $b}keys %sums;
        #@reads = sort { $a->{sum} <=> $b->{sum} } @reads;
        $self->warn_msg("Average base quality of extend read is lower than threshold 15, No extend")
            and return unless $avg{$order[$#order]} > 15;

		my $extend_read = $sums{pop @order};
			
		my $newbs = { type => 'base_segment', 
					  start => ($self->length + 1),
					  stop => $extend_read->align_stop,
					  read_name => $extend_read->name };
		my @bs = @{$self->base_segments};
		push @bs, $newbs;
		$self->base_segments(\@bs);
		$extend_read->align_clip_stop($extend_read->length);
		$extend_read->qual_clip_stop($extend_read->length);
		my $length = $extend_read->align_stop - $self->length;
		my $read_offset = $extend_read->length - $length;
		my $consensus = $self->padded_base_string.substr($extend_read->padded_base_string,$read_offset,$length);
		my @temp_quality = @{$extend_read->padded_base_quality};
		@temp_quality = splice(@temp_quality,$read_offset,$length);
		unshift @temp_quality,@{$self->padded_base_quality};
		$self->padded_base_string($consensus);
		$self->padded_base_quality(\@temp_quality);			
	}
    return 1;
}


sub get_phd {
    my ($phd_name, $po) = @_;
    my $err;

    for (@$po) {
        my $phd;
        eval {$phd = $_->get_phd($phd_name)};
        $err = $@;
        return $phd if $phd and !$err;
    }

    die $err if $err;
    die "can not find phd for $phd_name\n";
}
            
    

1;

#$HeadURL$
#$Id$
