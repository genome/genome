package Finishing::Assembly::ContigTools;
our $VERSION = 0.01;

use strict;

use Finishing::Assembly::Transform;
use Carp::Assert;
use Cwd;
use Bio::Seq;
use Bio::Tools::dpAlign;
use Finishing::Assembly::Utility;
use Gtk2;
use GSC::IO::Assembly::Ace::Writer;
#use Finishing::Assembly::Ace::Reader;
#use Finishing::Assembly::Ace;
use List::Util qw(min max);
use Finishing::Assembly::Factory;
use Finishing::Assembly::PhdDB;
use Finishing::Assembly::Phd;
#use GSC::Sequence::Assembly::AceAdaptor;
use Finishing::Assembly::AlignBP;
use Finishing::Assembly::AlignBPUP;
use Finishing::Assembly::AlignDerive;
use Data::Dumper;
use Finfo::Logging;

my $DEBUG = 0;

sub new 
{
    croak("__PACKAGE__:new:no class given, quitting") if @_ < 1;
	my ($caller, %params) = @_; 
    my $caller_is_obj = ref($caller);
    my $class = $caller_is_obj || $caller;
    my $self = {};
    bless ($self, $class);    
    return $self;
}

sub _align_contigs
{
	my ($left_contig,$right_contig,$gap_open_penalty,$gap_ext_penalty) = @_;
	my $left_string = $left_contig->padded_base_string;
	$left_string =~ tr/*/r/;
	my $right_string = $right_contig->padded_base_string;
	$right_string =~ tr/*/r/;

	#lets grab the last 1000 bases of the left contig, and the first 1000 bases
	# of the right contig.  We only need to do this if they are longer than
	# 1000 bases, which is why I'm checking

	if(length($left_string) >= 1000)
	{
		$left_string = substr($left_string, length($left_string)-1000, 1000);	
	}

	if(length($right_string) >= 1000)
	{
		$right_string = substr($right_string, 0, 1000);	
	}

	my $left_seq = Bio::Seq->new( -display_id => $left_contig->name,
		-seq => $left_string

	);

	my $right_seq = Bio::Seq->new( -display_id => $right_contig->name,
		-seq => $right_string

	);

	my $seqstring = $right_seq->seq();
	$right_seq->alphabet('dna');
	$left_seq->alphabet('dna');
	my $factory = Bio::Tools::dpAlign->new ( -match => 1,
		-mismatch => -1,
		-gap => $gap_open_penalty,
		-ext => $gap_ext_penalty,
		-alg => Bio::Tools::dpAlign::DPALIGN_LOCAL_MILLER_MYERS);		

	my $out = $factory->pairwise_alignment($left_seq, $right_seq);
	_convert_bioperl_to_gsc_transform($out,$left_contig,$right_contig);
}

sub _convert_bioperl_to_gsc_transform
{
	my($out,$left_contig,$right_contig) = @_;

	#we'll go ahead and create a hash named 'align' which will contain all of the
	#data that we are going to use when performing our alignment
	#think of it as a convenient way of passing around the alignment data for each
	#contig
	$left_contig->align({});
	$left_contig->align->{start}= $out->{_start_end_lists}{$left_contig->name}[0]{start};
	$left_contig->align->{end} = $out->{_start_end_lists}{$left_contig->name}[0]{end};

	#appears in the left contig to the one produced by bioperl
	#we use it later when adding these pads to the reads	
	$left_contig->align->{local_align_seq} = $left_contig->align->{align_seq} = $out->{_start_end_lists}{$left_contig->name}[0]{seq};
	$right_contig->align->{local_align_seq} = $right_contig->align->{align_seq} = $out->{_start_end_lists}{$right_contig->name}[0]{seq};
	$left_contig->align->{align_seq} =~ tr/rR/**/;
	$right_contig->align->{align_seq} =~ tr/rR/**/;
	if( $left_contig->length > 1000)
	{
		$left_contig->align->{start} = $left_contig->length-1000+$left_contig->align->{start};
		$left_contig->align->{end} = $left_contig->length-1000+$left_contig->align->{end};
	}
	else
	{
		$left_contig->align->{start} =	$left_contig->align->{start};
		$left_contig->align->{end} = $left_contig->align->{end};
	}
	$right_contig->align->{start}= $out->{_start_end_lists}{$right_contig->name}[0]{start};
	$right_contig->align->{end} = $out->{_start_end_lists}{$right_contig->name}[0]{end};
	my $left_pre = substr( $left_contig->padded_base_string, 0, $left_contig->align->{start}-1 );
	my $left_post = substr( $left_contig->padded_base_string, $left_contig->align->{end}, $left_contig->length-$left_contig->align->{end});
	$left_contig->align->{merge_seq} = $left_pre.$left_contig->align->{align_seq};
	$left_contig->align->{align_seq} = $left_pre.$left_contig->align->{align_seq}.$left_post;	
	my $right_pre = substr( $right_contig->padded_base_string, 0, $right_contig->align->{start}-1 );
	my $right_post = substr( $right_contig->padded_base_string, $right_contig->align->{end}, $right_contig->length-$right_contig->align->{end});
	$right_contig->align->{align_seq} = $right_pre.$right_contig->align->{align_seq}.$right_post;
	$left_contig->align->{merge_seq} .= $right_post;
	$left_contig->align->{merge_seq} =~ tr/-/*/;


	my $left_transform = Finishing::Assembly::Transform->new($left_contig->align->{align_seq}, '-');#this creates a transform for the merge region as it originally 
	my $right_transform = Finishing::Assembly::Transform->new($right_contig->align->{align_seq}, '-');#this creates a transform for the merge region as it originally 
	$right_transform->_offset( $left_contig->align->{start}-$right_contig->align->{start});
	#may want to check if offset is valid
	#go ahead and convert pads
	$left_contig->align->{local_align_seq} =~ tr/Rr-/***/;
	$right_contig->align->{local_align_seq} =~ tr/Rr-/***/;
	$left_contig->align->{transform} = $left_transform;
	$right_contig->align->{transform} = $right_transform;
}

sub _clip_read_after_position
{
	my ($read, $transform,$position) = @_;
	my $clip_position = $read->get_child_position_from_parent_position($position);
	$read->align_clip_start(min(max($read->align_clip_start,$clip_position),$read->align_stop));
	$read->qual_clip_start(min(max($read->qual_clip_start,$clip_position),$read->align_stop));
	$read->qual_clip_stop(min(max($read->qual_clip_stop,$clip_position),$read->align_stop));
	$read->align_clip_stop(min(max($read->align_clip_stop,$clip_position),$read->align_stop));
		
	$read->align_clip_start($transform->get_pad_position($read->align_clip_start));		
	$read->align_clip_stop($transform->get_pad_position($read->align_clip_stop));		
	$read->qual_clip_start($transform->get_pad_position($read->qual_clip_start));		
	$read->qual_clip_stop($transform->get_pad_position($read->qual_clip_stop));
	if($read->align_clip_start>=$read->align_clip_stop||
	   $read->qual_clip_start>=$read->qual_clip_stop )
	{
		$read->align_clip_start($read->align_clip_stop($read->qual_clip_start($read->qual_clip_stop(1))));	
	}
}

sub _clip_read_before_position
{
	my ($read, $transform,$position) = @_;
	my $clip_position = $read->get_child_position_from_parent_position($position);
	$read->align_clip_start(max(min($read->align_clip_start,$clip_position)),1);
	$read->qual_clip_start(max(min($read->qual_clip_start,$clip_position)),1);
	$read->qual_clip_stop(max(min($read->qual_clip_stop,$clip_position)),1);
	$read->align_clip_stop(max(min($read->align_clip_stop,$clip_position)),1);
		
	$read->align_clip_start($transform->get_pad_position($read->align_clip_start));		
	$read->align_clip_stop($transform->get_pad_position($read->align_clip_stop));		
	$read->qual_clip_start($transform->get_pad_position($read->qual_clip_start));		
	$read->qual_clip_stop($transform->get_pad_position($read->qual_clip_stop));
	if($read->align_clip_start>=$read->align_clip_stop||
	   $read->qual_clip_start>=$read->qual_clip_stop)
	{
		$read->align_clip_start($read->align_clip_stop($read->qual_clip_start($read->qual_clip_stop(1))));	
	}

}

sub _create_merge_tag
{
	my ($contig_name,$left_end_of_merge, $right_end_of_merge) = @_;
	my @temptime = localtime;
	my $year = sprintf("%02d", $temptime[5] % 100);
	$year = "20$year";
	my $tempstring = sprintf("%4d-%02d-%02d %02d:%02d:%02d", $year, ($temptime[4]+1), $temptime[3], 
	$temptime[2], $temptime[1], $temptime[0] );
	my $factory = Finishing::Assembly::Factory->connect("source");
	return $factory->create_consensus_tag(parent => $contig_name, type => 'comment', start => $left_end_of_merge, stop => $right_end_of_merge, date => $tempstring,
	   source => "cmt", no_trans => "NoTrans");#, scope => "ACE");
}

sub _translate_tag_positions
{
	my ($right_contig,$merge_contig) = @_;#from_contig,to_contig
	my @tags = @{$right_contig->tags};
	#foreach (@tags){$_ = $_->copy;};#make sure we're working on a copy of the tags
	foreach my $tag (@tags)
	{
		# first convert the tag to new consensus coordinates, this will adjust for any pads that
		# are inserted into the merge region
		$tag->start ( $right_contig->align->{transform}->get_pad_position($tag->start, 0));
		$tag->stop ( $right_contig->align->{transform}->get_pad_position($tag->stop, 0));							
	}
	return @tags;
}

sub _set_tags_parent_name
{
	my ($parent_name, $tags) = @_;
	foreach my $tag (@{$tags})
	{
		$tag->parent($parent_name);	
	}
	return $tags;
}


sub compute_stats
{
	my ($left_contig, $right_contig) = @_;
	my %stats;
	#now we compute percent identity
	{
		$stats{start} = $left_contig->align->{transform}->get_pad_position($left_contig->align->{start});
		$stats{stop} = $left_contig->align->{transform}->get_pad_position($left_contig->align->{end});		
		$stats{left_contig} = $left_contig;
		$stats{right_contig} = $right_contig;
	}
	
	$stats{mismatch} =0;
	$stats{length} = length($left_contig->align->{local_align_seq});
	for(my $pos=0;$pos<$stats{length};$pos++)
	{
		if(lc(substr($left_contig->align->{local_align_seq},$pos,1)) ne lc(substr($right_contig->align->{local_align_seq}, $pos, 1)))
		{
			$stats{mismatch}++;
		}
	}
	$stats{percent_identity} = 0;
	if($stats{length}>0)
	{
		$stats{percent_identity} = 100.0 - ((100.0*$stats{mismatch})/$stats{length});
	}
	else
	{
		$stats{percent_identity} = 0;
	}

	#compute high quality percentage match
	my @qarray1;
	for(my $pos = $left_contig->align->{start};$pos <= $left_contig->align->{end};$pos++)
	{
		push @qarray1, $left_contig->padded_base_quality->[$pos-1];
	}
	my @qarray2;
	for(my $pos = $right_contig->align->{start};$pos <= $right_contig->align->{end};$pos++)
	{
		push @qarray2, $right_contig->padded_base_quality->[$pos-1];
	}
	$stats{hqcount} = 0;
	$stats{hqmismatch}=0;
	
	for(my $pos=0;$pos<$stats{length};$pos++)
	{
		if($qarray1[$pos]>25&&$qarray2[$pos]>25)
		{
			$stats{hqcount}++;
			if(lc(substr($left_contig->align->{local_align_seq},$pos,1)) ne lc(substr($right_contig->align->{local_align_seq}, $pos, 1)))
			{
				$stats{hqmismatch}++;
				$stats{match_string} .= "H";
			}
			else
			{
				$stats{match_string} .= " ";
			}
		}
		else
		{
			if(lc(substr($left_contig->align->{local_align_seq},$pos,1)) ne lc(substr($right_contig->align->{local_align_seq}, $pos, 1)))
			{
				$stats{match_string} .= "D";
			}
			else
			{
				$stats{match_string} .= " ";
			}
		}
	}


	if($stats{hqcount}>0)
	{
		$stats{hq_percent_identity} = 100.0 - ((100.0*$stats{hqmismatch})/$stats{hqcount});
	}
	else
	{
		$stats{hq_percent_identity} = 0;
	}
	return \%stats;
}

sub print_stats2
{
	my ($stats, $statsfh) = @_;
	my $fh = $statsfh;
	if(!defined $fh){$fh =\*STDOUT};
	#now we compute percent identity
	my $left_contig = $stats->{left_contig};
	my $right_contig = $stats->{right_contig};
		
    print $fh "Merging ",$left_contig->name, " and ", $right_contig->name, " from consensus position $stats->{start} to $stats->{stop}. \n";

    #print $fh "Merge sequence is:\nLeft Contig:  $left_contig->align->{local_align_seq}\n";
    #print $fh "Right Contig: $right_contig->align->{local_align_seq}\n";
		
    #print $fh "Discrepancy:  $stats->{match_string}\n\n";
	print $fh "Statistics Info:\n";
	print $fh "Length of merge region:               $stats->{length}\n";
	print $fh "Number of discrepancies:              $stats->{mismatch}\n";
	printf $fh ("Percent Identity:                     %2.2f\%\n", $stats->{percent_identity});
	print $fh "Number of high quality bases:         $stats->{hqcount}\n";
	print $fh "Number of high quality discrepancies: $stats->{hqmismatch}\n";

	printf $fh ("High Quality Base Percent Identity:   %2.2f\%\n\n", $stats->{hq_percent_identity});

}

sub check_stats
{
	my ($stats, $check) = @_;
	return 1 if (!defined $check||!defined $stats);
	if (defined $check->{hqlength}&&$stats->{hqcount}<$check->{hqlength})
	{
		print STDERR "alignment region hq length cutoff not met.\n"; return 0;
	}
	if (defined $check->{hqmismatch}&&$stats->{hqmismatch}<$check->{hqmismatch})
	{
		print STDERR "hqmismatch cutoff not met.\n"; return 0;
	}
	if(defined $check->{hq_percent_identity}&&$stats->{hq_percent_identity}<$check->{hq_percent_identity})
	{
		print STDERR "hq_percent_identity cutoff not met.\n"; return 0;
	}
	if (defined $check->{percent_identity}&&$stats->{percent_identity}<$check->{percent_identity})
	{
		print STDERR "percent_identity cutoff not met.\n"; return 0;
	}
	if (defined $check->{mismatch}&&$stats->{mismatch}<$check->{mismatch})
	{
		print STDERR "mismatch cutoff not met.\n"; return 0;
	}
	if (defined $check->{length}&&$stats->{length}<$check->{length})
	{
		print STDERR "alignment region length cutoff not met.\n"; return 0;
	}	
	return 1;		
}

sub print_stats
{
	my ($left_contig, $right_contig) = @_;

	#now we compute percent identity
	{
		my $temp1 = $left_contig->align->{transform}->get_pad_position($left_contig->align->{start});
		my $temp2 = $left_contig->align->{transform}->get_pad_position($left_contig->align->{end});
		print "Merging ",$left_contig->name, " and ", $right_contig->name, " from consensus position $temp1 to $temp2. \n";
	}
	print "Merge sequence is:\nLeft Contig:  $left_contig->align->{local_align_seq}\n";
	print "Right Contig: $right_contig->align->{local_align_seq}\n";
	my $mismatch =0;
	my $length = length($left_contig->align->{local_align_seq});
	for(my $pos=0;$pos<$length;$pos++)
	{
		if(lc(substr($left_contig->align->{local_align_seq},$pos,1)) ne lc(substr($right_contig->align->{local_align_seq}, $pos, 1)))
		{
			$mismatch++;
		}
	}
	my $percent_identity;
	if($length>0)
	{
		$percent_identity = 100.0 - ((100.0*$mismatch)/$length);
	}
	else
	{
		$percent_identity = 0;
	}

	#compute high quality percentage match
	my @qarray1;
	for(my $pos = $left_contig->align->{start};$pos <= $left_contig->align->{end};$pos++)
	{
		push @qarray1, $left_contig->padded_base_quality->[$pos-1];
	}
	my @qarray2;
	for(my $pos = $right_contig->align->{start};$pos <= $right_contig->align->{end};$pos++)
	{
		push @qarray2, $right_contig->padded_base_quality->[$pos-1];
	}
	my $hqcount = 0;
	my $hqmismatch=0;
	my $sMatchString;
	for(my $pos=0;$pos<$length;$pos++)
	{
		if($qarray1[$pos]>25&&$qarray2[$pos]>25)
		{
			$hqcount++;
			if(lc(substr($left_contig->align->{local_align_seq},$pos,1)) ne lc(substr($right_contig->align->{local_align_seq}, $pos, 1)))
			{
				$hqmismatch++;
				$sMatchString .= "H";
			}
			else
			{
				$sMatchString .= " ";
			}
		}
		else
		{
			if(lc(substr($left_contig->align->{local_align_seq},$pos,1)) ne lc(substr($right_contig->align->{local_align_seq}, $pos, 1)))
			{
				$sMatchString .= "D";
			}
			else
			{
				$sMatchString .= " ";
			}
		}
	}

	my $nHqPercentIdentity;
	if($hqcount>0)
	{
		$nHqPercentIdentity = 100.0 - ((100.0*$hqmismatch)/$hqcount);
	}
	else
	{
		$nHqPercentIdentity = 0;
	}
	print "Discrepancy:  $sMatchString\n\n";
	print "Statistics Info:\n";
	my $templength = length( $left_contig->align->{local_align_seq} );
	print "Length of merge region:               $templength\n";
	print "Number of discrepancies:              $mismatch\n";
	printf ("Percent Identity:                     %2.2f\%\n", $percent_identity);
	print "Number of high quality bases:         $hqcount\n";
	print "Number of high quality discrepancies: $hqmismatch\n";

	printf ("High Quality Base Percent Identity:   %2.2f\%\n", $nHqPercentIdentity);
}

sub get_phd
{
	my ($phd_name, $output_phd, @phd_array) = @_;
	my $error;
    my $bac_phd = 1;
    if ($phd_name =~ /aae13h12/){
        my $var = 1;
    }
    foreach (@phd_array)
    {
		my $phd;
		eval
		{
			$phd = $_->get_phd($phd_name);
            if(0)
            #if(defined $output_phd)
            {
                $output_phd->add_phd($phd) if $phd and $bac_phd ==1;
            }
            $bac_phd =0;
        };
        $error = $@;
		return $phd if $phd and !$error;
		
	}
	if ($error) {die "$error"};
    die "Could not find phd";
}

sub merge
{
	my ($self, $left_contig, $right_contig, $phd_object, %params) = @_;
	
	my $output_phd = $params{output_phd};
	my $nGapOpenPenalty=$params{gap_open_penalty} || 15;
	my $nGapExtPenalty=$params{gap_ext_penalty} || 15;
	my $nGGapExtPenalty=$params{glob_gap_open_penalty} || 15;
	my $nGGapOpenPenalty=$params{glob_gap_open_penalty} || 15;
	my $bUseGlobalAlign=$params{use_global_align} || 0;
	my $bQuiet=$params{quiet} || 0;	
	my $phd_array=$params{phd_array};
	my @phd_array = ($phd_object) if(defined $phd_object);
	push @phd_array,@{$phd_array} if(defined $phd_array);
	
	my %left_reads = map {$_->name,$_ } $left_contig->assembled_reads->all;
	my %right_reads = map {$_->name,$_ } $right_contig->assembled_reads->all;
	my %left_overlap_reads;
	my %right_overlap_reads;
	my @base_segs;		
    $self->info_msg("1");
	assert($left_contig->is_bs_array_structure_ok(1));
    $self->info_msg("2");
	assert($right_contig->is_bs_array_structure_ok(1));
    $self->info_msg("3");

	#_align_contigs($left_contig,$right_contig,3,5);
	my $max_length;
	$max_length = 3000 if (!exists($params{max_length}));
	%params = (%params,gap_open_penalty => 4, gap_ext_penalty => 2, max_length => $max_length);
	my @align;
	if(!exists $params{derive_read})
	{
		@align = Finishing::Assembly::AlignBP::align($left_contig->padded_base_string, $right_contig->padded_base_string, %params);	
	}
	else
	{
		$params{read_name} = $params{derive_read};#'fake2.b1';#'TPAC-30H12a.b1';#'fake2.b1';
		@align = Finishing::Assembly::AlignDerive::align($left_contig, $right_contig, %params);	
	}
    $self->info_msg("4");
	#my $dumper =  Dumper @align;
	#my $dumpfh = IO::File->new(">dumpeddata");
	#print $dumpfh $dumper;
	#exit;
	$left_contig->align( $align[0]);
	$right_contig->align( $align[1]);
	#print_stats($left_contig, $right_contig);
	my $stats = compute_stats($left_contig, $right_contig);
	print_stats2($stats, $params{statsfh});
	if(!check_stats($stats,$params{cutoffs})){ return };	
	
	#take the left contig up to the last 1000 bases and repad it.

	#here is where we merge the consensus for the two contigs, basically, we know that the locally aligned region
	#will be the same for the left and right contigs.  So, we just add the consensus in the right contig that occurs
	#after the locally aligned region to the left contig.

	#my $merge_contig = Finishing::Assembly::Contig->new;
	my $factory = Finishing::Assembly::Factory->connect("source");
	my $merge_contig = $factory->create_contig(name => $left_contig->name);
    $self->info_msg("5");
	$merge_contig->name( $left_contig->name );
	$merge_contig->align( {});
	$merge_contig->align->{start} = $left_contig->align->{transform}->get_pad_position($left_contig->align->{start});
	$merge_contig->align->{end} = $left_contig->align->{transform}->get_pad_position($left_contig->align->{end});
	$merge_contig->complemented( $left_contig->complemented );
	$merge_contig->padded_base_string($left_contig->align->{merge_seq});	
    $self->info_msg("6");
	
	foreach(values %left_reads)
	{
		
		if($_->align_stop > $left_contig->align->{start} )
		{	
			my $new_seq = $left_contig->align->{transform}->pad_string_partial($_->padded_base_string, $_->position-1, '-');
			my $transform = Finishing::Assembly::Transform->new($new_seq, '-');
			$new_seq =~ tr/-/*/;
			$_->padded_base_string( $new_seq);


			# set clipping to before the end of the merge region if necessary
			_clip_read_before_position($_,$transform,$left_contig->align->{end});
			#$_->{padded_length} = length($_->padded_base_string);
			#$_->{fromcontig} = "left";
			$left_overlap_reads{$_->name()} = $_;
		}
		$_->position ($left_contig->align->{transform}->get_pad_position($_->position));
	}
    $self->info_msg("7");
	delete $right_reads{$params{derive_read}} if (exists $params{derive_read});
	foreach(values %right_reads)
	{
		if($_->position <= $right_contig->align->{end} )
		{
			my $tempoffset = 0;#$right_contig->align->{transform}->_offset;
			my $new_seq = $right_contig->align->{transform}->pad_string_partial($_->padded_base_string, $_->position-$tempoffset-1, '-');
			my $transform = Finishing::Assembly::Transform->new($new_seq, '-');
			$new_seq =~ tr/-/*/;
			$_->padded_base_string( $new_seq);


			# set clipping to before the end of the merge region if necessary
			_clip_read_after_position($_,$transform,$right_contig->align->{start});

			#$_->{padded_length} = length($_->padded_base_string);
			#$_->{fromcontig} = "right";
			$right_overlap_reads{$_->name()} = $_;
		}
		$_->position ($right_contig->align->{transform}->get_pad_position($_->position));
	}
    $self->info_msg("8");

	#go ahead and add the reads together
	#add the left reads to the list of merge reads
	my @overlap_reads;
	foreach( values %left_overlap_reads,values %right_overlap_reads)
	{

		my $phd = get_phd($_->phd_file, $output_phd, @phd_array);
		my @quality = @{$phd->unpadded_base_quality};
		@quality = reverse @quality if ($_->complemented);
		$_->unpadded_base_quality( \@quality);
		push @overlap_reads, $_;		
	}
	$merge_contig->assembled_reads([values %left_reads, values %right_reads]);
	
	#now that we have our list of reads in the merge region, lets open the quality values
	#for these reads so that we can recalculate the quality values for the consensus

	my $left_end_of_merge = $merge_contig->align->{start};
	my $right_end_of_merge = $merge_contig->align->{end};
		
	# reset the consensus base and quality arrays to refigure them
	# base on the combined reads from both the left and right contig.  Note that it is
	# possible (likely) that the consensus bases and quality in the
	# merged region will be different for the new contig.  
	
	my ($best_quality_reads,$best_quality) = $merge_contig->recalculate_consensus($merge_contig->align->{start},
										 $merge_contig->align->{end},
										 \@overlap_reads);
	
    $self->info_msg("9");
	# append base segments for the left contig before the merge region
	my $bs_position = $left_contig->get_bs_index_at_position( $left_contig->align->{start});
	push( @base_segs, @{$left_contig->base_segments} );
	#don't add BS where there aren't any, make sure the array is at least one unit long
	if($base_segs[$bs_position]{start}<$merge_contig->align->{start})
	{
		#in this case, the base segment at this position starts before the merge region
		#clip the end pos to the base before the merge region
		$base_segs[$bs_position]{stop} = $merge_contig->align->{start}-1;	
		$#base_segs = $bs_position;
	}
	else
	{
		$#base_segs = $bs_position-1;
	}
	#assert( $merge_contig->is_bs_array_ok(0) );			

    $self->info_msg("10");
	my @bs_temp = @{$merge_contig->calculate_base_segments($merge_contig->align->{start},$merge_contig->align->{end},$best_quality_reads)};
	push(@base_segs, @bs_temp);
	
	#assert( $merge_contig->get_bs_array_structure_ok(0) );			
	# get the Base Segments that start after the merge region, and recalculate their start and end positions,
	# then push them onto the base seg array (after the merge region)
	$bs_position = $right_contig->get_bs_index_at_position( $right_contig->align->{end}+1);
	@bs_temp = @{$right_contig->base_segments};
	@bs_temp = splice ( @bs_temp, $bs_position, @bs_temp - $bs_position);
    $self->info_msg("11");
	if(@bs_temp>0&&$bs_temp[0]->{stop}>=($right_contig->align->{end}+1))
	{
		foreach my $bs (@bs_temp)
		{
			$bs->{start} = $right_contig->align->{transform}->get_pad_position($bs->{start});
			$bs->{stop} = $right_contig->align->{transform}->get_pad_position($bs->{stop});
		}
		#clip the first BS of aBSTemp so that it starts AFTER the merge
		#region!
		if($bs_temp[0]{start}<=$base_segs[@base_segs-1]{stop})
		{
			$bs_temp[0]{start} = $base_segs[@base_segs-1]{stop}+1;
		}
	
		push(@base_segs, @bs_temp);
		#$merge_contig->base_segments(\@base_segs);
		#assert( $merge_contig->is_bs_array_structure_ok(1)) ;
	
	}	
	$merge_contig->base_segments(\@base_segs);
    $self->info_msg("12");
	assert( $merge_contig->is_bs_array_structure_ok(1)) ;
    $self->info_msg("13");
	
	# next, splice in the newly computed quality values above to the left contig, then we will append the quality
	
	my @cons_quality = @{$left_contig->padded_base_quality};
	$#cons_quality = $left_contig->align->{start}-2;

	my @temp = @{$best_quality};
	splice(@temp,0,$merge_contig->align->{start});
	push (@cons_quality, @temp);

	my @right_cons_quality = @{$right_contig->padded_base_quality};	
	splice(@right_cons_quality,0,$right_contig->align->{end});
	push (@cons_quality,@right_cons_quality);
	$merge_contig->padded_base_quality(\@cons_quality);
    $self->info_msg("14");
	# now recalculate consensus qualities

	Finishing::Assembly::Utility::recalculate_consensus_qualities_and_change( $merge_contig, $merge_contig->align->{start}, $merge_contig->align->{end}, \@overlap_reads  );

	#transfer left contig tags and translated right contig tags
	my @merge_contig_tags = (@{$left_contig->tags}, _translate_tag_positions($right_contig,$merge_contig));
    $self->info_msg("15");
	push @merge_contig_tags, _create_merge_tag($merge_contig->name,$merge_contig->align->{start},$merge_contig->align->{end});	
    $self->info_msg("16");
	#_set_tags_parent_name($merge_contig->name, \@merge_contig_tags);
	$merge_contig->tags(\@merge_contig_tags	);
    $self->info_msg("17");
	
	return $merge_contig;
}

	


sub _get_ace_object
{
	my ($self, $data_source, @contig_names) = @_;
	my $in_fh;
	#if(!($data_source =~ /\.ace/))
#	{
#		my $assembly;# = GSC::Sequence::Assembly->get( sequence_item_name => $data_source);
#
#		foreach (@contig_names)
#		{
#			$_ = "$data_source.$_";
#		}
#		my @ctgs = GSC::Sequence::Contig->get( sequence_item_name => [@contig_names]);
#		foreach my $ctg (@ctgs) {
#		   $ctg->lock;
#		   #print $ctg->sequence_item_name,"\n";
#		}
#		my $contig_names = join ' ', @contig_names;
#
#		print "Grabbing $contig_names from the database\n";
#
#		my $in_fs;
#		$in_fh = new IO::String($in_fs);
#		my $ace_writer = GSC::IO::Assembly::Ace::Writer->new($in_fh);
#		my $ace_adapter = GSC::Sequence::Assembly::AceAdaptor->new();
#
#		$ace_adapter->export_assembly(
#			assembly => $assembly,
#			writer => $ace_writer,
#			contigs => \@ctgs
#		);
#
#		$in_fh->setpos(0);	
#		
#	}
#	else
#	{
#		$in_fh = new IO::File($data_source);
#		
#	}
	#return Finishing::Assembly::Ace->new(input => $in_fh);
}

sub _get_phd_object
{
	my ($self, $data_source) = @_;
	if(!($data_source =~ /\.ace/))
	{		
		return Finishing::Assembly::PhdDB->new;
	}
	else		
	{
		my $cwd = getcwd;
		return Finishing::Assembly::Phd->new(input_directory => "$cwd/../phd_dir/");		
	}
}

sub _calculate_split_region
{
	my ($left_reads, $right_reads, $contig_length) = @_;
	my $rRightMostLeftRead;
	my $rLeftMostLeftRead;
	my $nLeftMostLeftRead;
	my $nRightMostLeftRead;

	my $init = 0;
	foreach my $read (@{$left_reads})
	{
		if ( $init == 0 )
		{
			$nRightMostLeftRead = $read->align_stop;
			$rRightMostLeftRead = $read;
			$nLeftMostLeftRead  = $read->align_start;
			$rLeftMostLeftRead = $read;
			$init = 1;
		}
		else
		{
			if ( $nRightMostLeftRead < $read->align_stop )
			{
				$nRightMostLeftRead = $read->align_stop;
				$rRightMostLeftRead = $read;
			}
			if ( $read->align_start < $nLeftMostLeftRead )
			{
				$nLeftMostLeftRead = $read->align_start;
				$rLeftMostLeftRead = $read;
			}
		}
	}

	my $rLeftMostRightRead;
	my $rRightMostRightRead;
	my $nRightMostRightRead;
	my $nLeftMostRightRead;
	$init = 0;

	foreach my $read (@{$right_reads})
	{
		if ( $init == 0 )
		{
			$nLeftMostRightRead = $read->align_start;
			$rLeftMostRightRead  = $read;
			$nRightMostRightRead = $read->align_stop;
			$rRightMostRightRead = $read;
			$init = 1;
		}
		else
		{
			if ( $read->align_start < $nLeftMostRightRead )
			{
				$nLeftMostRightRead = $read->align_start;
				$rLeftMostRightRead = $read;
			}
			if ( $nRightMostRightRead < $read->align_stop )
			{
				$nRightMostRightRead = $read->align_stop;
				$rRightMostRightRead = $read;
			}
		}
	}

	# now find the overlap region
	#
	# the assertion is due to the fact that, since the assembly was
	# original contiguous, there is no way that it can be divided into
	# two groups of reads that don't overlap

	# the start_, nRightEndOfTear_ region marks the torn region--the
	# region in which each column has a different number of reads than
	# it did in the original contig.  Note that this is this region:

	#	------
	#	   ---------------------------
	#		 ---------------------------
	#			 +++++++++++++++++++++++++++++++++++
	#				   +++++++++++++++++++++++++++++++++++
	#									  +++++++++++++++++++++++++++++++++++
	#			 ^  				   ^
	#			 |  				   |
	#			 start_	   nRightEndOfTear_
	# It is NOT:
	#		^											 ^
	#		here			  to						 here
	# That is, a read that is at the cursor where the user is tearing,
	# may not be completely contained within the torn region.
	
	my ($start, $stop);
	assert(
			calculate_intersection(
						$nLeftMostRightRead, $nRightMostRightRead,
						$nLeftMostLeftRead,  $nRightMostLeftRead,
						\$start,    \$stop
			)
	);
	
	# check that this is not happening:
	#  -------------------------------------
	#		  ++++++++++++++++
	# where --- goes into the left contig and +++ goes into the right
	# This would be a problem when moving base segments, so we just
	# disallow it and make the user use PutReadIntoItsOwnContig

	if ( $nLeftMostLeftRead > $nLeftMostRightRead )
	{
		printf("There is a problem because the new right contig has a read %s which has left position %d which is more left than any read in the new left contig.  (The leftmost position of the left contig is %d due to read %s.) Thus this is not a tear.  Instead you could use PutReadIntoItsOwnContig",
			$rLeftMostRightRead->name, $nLeftMostRightRead,
			$nLeftMostLeftRead,          $rLeftMostLeftRead->name
		);
		exit;		
	}

	if ( $nRightMostRightRead < $nRightMostLeftRead )
	{
		printf("There is a problem because the new left contig has a read %s which protrudes to the right to %d, further than any read in the new right contig.  (The rightmost read in the right contig is %s which extends to %d.)  Thus this is not a tear.  Instead you could use PutReadIntoItsOwnContig",
			$rRightMostLeftRead->name,  $nRightMostLeftRead,
			$rRightMostRightRead->name, $nRightMostRightRead
		);
		exit;		
	}
	
	#clip to within the contig region, bug fix
	$start = 1 if($start < 1);
	$stop = $contig_length if($stop > $contig_length);
	return ($start, $stop);
}

sub get_user_selected_reads
{
	my ($ref_reads) = @_;
	Gtk2->set_locale; 
	Gtk2->init;

	my %reads = %{$ref_reads};

	my $d = Gtk2::Dialog->new();
	$d->set_usize(300,300);
	$d->show;

	my $sw = Gtk2::ScrolledWindow->new();
	$sw->set_policy("automatic", "automatic");
	#$sw->set_policy("always", "always");
	$sw->show;
	$d->vbox->pack_start_defaults($sw);

	my $vbox = Gtk2::VBox->new();
	$vbox->show;
	$sw->add_with_viewport($vbox);

	foreach my $read_name (keys %reads)
	{
    	my $bc = Gtk2::ButtonCrate::Radio->new('h');
    	$vbox->pack_start($bc->bbox, 0, 0, 0);

    	my $label = Gtk2::Label->new($read_name);
    	$label->show;
    	$bc->bbox->pack_start($label, 0, 0, 0);

    	$bc->add_buttons(" <-", " ->");

		$bc->set_button_active_by_num(1) if($reads{$read_name} eq "right");	

    	$bc->bbox->set_usize(100, 50);

    	$reads{$read_name} = $bc;
	}

	my $b = Gtk2::Button->new("Done");
	$b->signal_connect
	(
    	"clicked",
    	sub
    	{
            my ($self, $reads, $ref_reads, $d) = @_;
            foreach my $read (keys %{$reads})
        	{
                my $bc = ${$reads}{$read};

            	#print $bc->get_active_button_name,"\n";
				if($bc->get_active_button_name eq " <-")
				{
					${$ref_reads}{$read} = "left";
				}
				else
				{
					${$ref_reads}{$read} = "right";
				}				
        	}
			$d->hide;
        	Gtk2->main_quit;
    	},
		\%reads,
		$ref_reads,
		$d
	);
	$b->show;
	$d->action_area->pack_start_defaults($b);

	$b = Gtk2::Button->new("Cancel");
	$b->signal_connect("clicked", sub{ Gtk2->exit(0) });
	$b->show;
	$d->action_area->pack_start_defaults($b);

	Gtk2->main;

}
sub split
{
	my ($self, $contig, $phd_object, %params) = @_;		

	my $split_position = $params{split_position} if exists $params{split_position};
	my $phd_array=$params{phd_array};
	my @phd_array = ($phd_object) if(defined $phd_object);
	push @phd_array,@{$phd_array} if(defined $phd_array);	
	my $in_fh;
	my $out_fh;

	my $bKeepCounting=1;
	my $no_gui=0;
	$no_gui = $params{'no_gui'} if exists $params{'no_gui'};		

	my %user_select_reads;
    #my $padded_split_position = $split_position;#$contig->_transform->get_pad_position($split_position);
    my $padded_split_position = $contig->unpad_position_to_pad_position($split_position);

	my $reads = $contig->assembled_reads;	
	my @left_reads;
	my @right_reads;
	foreach my $read ($reads->all) 
	{
		last if (defined $params{reads_given});
		#	check whether this read is in the left or the right contig
		#   If it doesn't intersect the cursor location, but is on its left,
		#   it is in the left contig(s)

		#   If it doesn't intersect the cursor location, but is on its right,
		#   it is in the right contigs(s).

 		#   If it *does* intersect the cursor location, then go by the
 		#   user's selection.  First we need to find the Read's position in the
 		#   consensus, which is why we convert from Read Index to Consensus coordinates

		my $read_position = ($read->get_parent_position_from_child_position ($read->qual_clip_start, 1) + $read->get_parent_position_from_child_position($read->qual_clip_stop, 1))/2.0;
		if($no_gui)
		{
			if ( $read_position < $padded_split_position )
			{
				push @left_reads, $read;
			}
			else
			{
				push @right_reads, $read;
			}
		}
		else
		{
			if ( $read->align_stop < $padded_split_position)
			{
				push @left_reads, $read;
			}
			elsif( $read->align_start > $padded_split_position )
			{
				push @right_reads, $read;
			}
			else
			{
				if($read_position < $padded_split_position)
				{
					$user_select_reads{$read->name} = "left";
				}
				else
				{
					$user_select_reads{$read->name} = "right";
				}	
			}
		}
	}
	if(defined $params{reads_given})
	{
		@left_reads = @{$params{left_reads}};
		@right_reads = @{$params{right_reads}};
	}

	if(keys %user_select_reads)
	{
		get_user_selected_reads(\%user_select_reads);
		my %reads = map {$_->name,$_} $reads->all;
		foreach my $read_name (keys %user_select_reads)
		{
			if($user_select_reads{$read_name} eq "left")
			{
				push @left_reads, $reads{$read_name};
			}
			else
			{
				push @right_reads, $reads{$read_name};
			}	
		}
	}
	assert( ( scalar @left_reads + scalar @right_reads ) == $reads->count );
	if($params{allow_no_split} == 1)
	{
		if( (scalar @left_reads) == 0)
		{
			return (undef, $contig);
		}
		if( (scalar @right_reads) == 0)
		{
			return ($contig, undef);
		}
	}
	else
	{	
		die(   "You are putting all of the reads into the right contig so this is not a tear") if ( (scalar @left_reads) == 0 ) ;
		die(   "You are putting all of the reads into the left contig so this is not a tear") if ( (scalar @right_reads) == 0 );
	}
	# at this point, all reads in the contig have been put into 2 lists:
	# a list of reads that the user wants to go into the 'left' contig,
	# and a list of reads that the user wants to go into the 'right contig.
	# (Note that these contigs themselves may fall apart into more contigs.)

	# Find the torn region--To the left and right of this region, any column
	# will be precisely the same as it was before the tear:  the same reads
	#  be in that column.  Within the tear, the reads will be different,
	# and thus we will have to recompute the consensus, Golden path, and
	# base segment.

	my ($split_region_start, $split_region_stop) = _calculate_split_region(\@left_reads,\@right_reads,$contig->length);
    my $transform = Finishing::Assembly::Transform->new($contig->padded_base_string);
	print "tearing contigs from ", $transform->get_unpad_position($split_region_start), 
		  " to ", $transform->get_unpad_position($split_region_stop), "\n";

	
	# now lets find if things fall into more than 2 contigs in the torn region
	# To do this, let's make two little lists:  The left reads in the
	# torn region, and the right reads in the torn region.
	my $left_overlap_reads = $contig->get_overlapping_reads($split_region_start, $split_region_stop, \@left_reads);
	my $right_overlap_reads = $contig->get_overlapping_reads($split_region_start, $split_region_stop, \@right_reads);

	foreach( @$left_overlap_reads, @$right_overlap_reads)
	{
		my $phd = get_phd($_->phd_file, undef, @phd_array);
		my @quality = @{$phd->unpadded_base_quality};
		@quality = reverse @quality if ($_->complemented);
		$_->unpadded_base_quality( \@quality);				
	}
	# Now let's check (for the left contig and the right contig)
	# that there is a base for every position
	# (unaligned or aligned) for this entire region

	exit if ( !check_contiguous( $left_overlap_reads, $split_region_start, $split_region_stop ) || 
		 	  !check_contiguous( $right_overlap_reads, $split_region_start, $split_region_stop ));	
	
	# point of no return
	# At this point, there is no longer any error checking.
	# Errors indicate a program bug.

	# Transfer reads to new contigs.
	# Fix alignment positions and tag positions of the right contig 326-592, 592-870
	my $factory = Finishing::Assembly::Factory->connect("source");
	my $left_contig = $factory->create_contig(name => $contig->name . "a");
	#my $left_contig = Finishing::Assembly::Contig->new;
	#$left_contig->name( $contig->name . "a" );
	$left_contig->assembled_reads(\@left_reads);
	$left_contig->complemented( $contig->complemented);

	# in the torn region, figure out the consensus base and quality
	# The fast method of doing this is to go through each read once,
	# and mark two arrays--the highest quality at a cons pos and the
	# read that attained that highest quality.  I will ignore the
	# "aligned" clipping at this point, since a read could be high
	# quality unaligned precisely because there is a mis-join, which
	# is the reason the user is tearing the contig.  Thus the correct
	# consensus could be this unaligned high quality region.
	
	my $left_cont_start = max( $split_region_start,  1 );
	my $left_cont_stop   = min( $split_region_stop, $contig->length );
	assert( $left_cont_start <= $left_cont_stop );
	$left_contig->padded_base_string(substr( $contig->padded_base_string, 0, $left_cont_stop ));
	my @temp_qual = @{$contig->padded_base_quality};
	$left_contig->padded_base_quality([ splice(@temp_qual, 0, $left_cont_stop)]);
	my ($best_quality_reads, $best_quality) = $left_contig->recalculate_consensus($split_region_start,$split_region_stop,$left_overlap_reads);
		
	# transferring base segments

	# which is the last base segment that can be transferred?

	# There could be a situation in which we are tearing a contig like this:
	#			cccccccccccccccccccccccccc
	#	---------------------------
	#		 +++++++++++++++++++++++++++++++++++
	#
	# where the --- go into the new left contig and the
	# +++ go into the new right contig.
	# In this case, the torn region starts before the consensus
	# and thus there is no base segments before the torn region.
	# Thus no transferring of base segments before the torn region.
	assert( $contig->is_bs_array_structure_ok( 0 ));

	my @left_bs;     # Array of Left Base Segments
	my @right_bs;    # Array of Right Base Segments
	my @bs = @{$contig->base_segments};
	if ( 1 < $split_region_start )
	{
		my $nAfterLastBaseSegmentToBeTransferred =
		  $contig->get_bs_index_at_position( $split_region_start );

		if ( $nAfterLastBaseSegmentToBeTransferred == -1 )
		{
			print "nAfterLastBaseSegmentToBeTransferred == -1, start = ",
			  $split_region_start, "\n";
			assert(0);
		}

		my $nLastBaseSegmentToBeTransferred =
		  $nAfterLastBaseSegmentToBeTransferred - 1;

		# note that if the base segment that covers start_ is the
		# first base segment, then nLastBaseSegmentToBeTransferred will
		# be -1, so no base segments will be transferred in the loop
		# below.

		for ( my $nBaseSeg = 0 ;
			  $nBaseSeg <= $nLastBaseSegmentToBeTransferred ;
			  $nBaseSeg++ )
		{
			my $rBaseSeg = $bs[$nBaseSeg];

			my $rNewBaseSeg = {
								type      => 'base_segment',
								start => $rBaseSeg->{start},
								stop   => $rBaseSeg->{stop},
								read_name => $rBaseSeg->{read_name},								
			};
			push( @left_bs, $rNewBaseSeg );
		}

		# now look at the base segment that perhaps must be cut short
		# We now guarantee that this base segment (or the one we just added),
		# will end precisely before the tear region

		my $rBoundaryBaseSeg = $bs[ $nLastBaseSegmentToBeTransferred + 1 ];

		if ( $rBoundaryBaseSeg->{start} < $split_region_start )
		{
			assert( $rBoundaryBaseSeg->{stop} >= $split_region_start );
			my $rNewBoundaryBaseSeg = {
										type      => 'base_segment',
										start => $rBoundaryBaseSeg->{start},
										stop   => ( $split_region_start - 1 ),
										read_name => $rBoundaryBaseSeg->{read_name},
			};

			push( @left_bs, $rNewBoundaryBaseSeg );
		}
	}

	# now add base segments for the tear region
	# We will add a new base segment for each base and not worry about
	# trying to consolidate adjacent base segments that share the same
	# read.

	# Fixed Bug:  DG, October 4, 1000.  It is possible that this torn region
	# can extend beyond the end of the consensus.  Clearly it doesn't
	# make sense to have base segments there, in either new contig.  So
	# only create base segments in the torn region that overlaps the consensus.

	my $nTornRegionOnConsensusLeft;
	my $nTornRegionOnConsensusRight;

	if (
		calculate_intersection(
					$split_region_start,
					$split_region_stop,
					1,                              #start of consensus
					$left_contig->length,         #end of consensus
					\$nTornRegionOnConsensusLeft,
					\$nTornRegionOnConsensusRight
		)
	  )
	{
		for ( my $pos = $nTornRegionOnConsensusLeft ;
			  $pos <= $nTornRegionOnConsensusRight ;
			  $pos++ )
		{
			my $rNewBaseSegment = {
								   type      => 'base_segment',
								   start => $pos,
								   stop   => $pos,
								   read_name => $best_quality_reads->[$pos]->name,
			};
			push( @left_bs, $rNewBaseSegment );
		}
	}
	$left_contig->base_segments( \@left_bs);
	
	# now recalculate the consensus quality values for the new left contig
	# (The reason for the intersect is that there may be reads that
	# stick out on the left to negative consensus positions and some
	# of these may be in the right contig and thus the start_
	# may be to the left of where the new left contig starts.)

	my $nConsPosLeftToRecalculate;
	my $nConsPosRightToRecalculate;

	if (
		 calculate_intersection(
					 1,                      $left_contig->length,
					 $split_region_start,             $split_region_stop,
					 \$nConsPosLeftToRecalculate, \$nConsPosRightToRecalculate
		 )
	  )
	{
		Finishing::Assembly::Utility::recalculate_consensus_qualities_and_change($left_contig, $nConsPosLeftToRecalculate, $nConsPosRightToRecalculate, $left_overlap_reads);
	}

	# That completes the base segments for the new left contig

	# Now work on the new right contig.
	
	# add reads to new contig

	my $offset = -$split_region_start + 1;

	# handle case like this:
	#			cccccccccccccccccccccccccc
	#	---------------------------
	#		 +++++++++++++++++++++++++++++++++++

	# where the --- go into the new left contig and the
	# +++ go into the new right contig.
	# In this case, the beginning of the ccc... is 1 in
	# both the new left contig and the new right contig
	# Thus there is no translation.

	if ( $split_region_start < 1 )
	{
		$offset = 0;
	}

	#now we recalculate the coordinates for the reads that are going into the right
	#contig
	foreach my $read (@right_reads)
	{
		$read->position($read->position + $offset);
	}
	
	#my $right_contig = Finishing::Assembly::Contig->new;
	my $right_contig = $factory->create_contig(name => $contig->name . "a");
	#$right_contig->name( $contig->name . "b" );
	$right_contig->assembled_reads( \@right_reads);
	$right_contig->complemented( $contig->complemented);

	# reset the consensus base and quality arrays to refigure them
	# based on the reads from the right contig.  Note that it is
	# possible (likely!) that the consensus bases and quality in the
	# torn region will be different for the two ends of the new contigs

	my $right_cont_start = max( $split_region_start,  1 );
	my $right_cont_stop = min( $split_region_stop, $contig->length );
	assert( $right_cont_start <= $right_cont_stop );
	$right_contig->padded_base_string(substr( $contig->padded_base_string, $right_cont_start-1, $contig->length-$right_cont_start + 1 ));
	@temp_qual = @{$contig->padded_base_quality};
	$right_contig->padded_base_quality([ splice(@temp_qual, $right_cont_start-1, $contig->length-$right_cont_start + 1)]);
	
	$right_cont_start += $offset;
	$right_cont_stop += $offset;
	($best_quality_reads, $best_quality) = $right_contig->recalculate_consensus($right_cont_start,$right_cont_stop,$right_overlap_reads);

	# now create base segments for the torn region

	# Fixed Bug:  DG, October 4, 1000.  It is possible that this torn region
	# can extend beyond the end of the consensus.  Clearly it doesn't
	# make sense to have base segments there, in either new contig.  So
	# only create base segments in the torn region that overlaps the consensus.
	# I saw this in a case in which the torn region went all the way to
	# the left end of the contig, before the consensus:
	#		   cccccccccccccccc   (consensus)
	#	lllllllllllllllll	(goes into left contig)
	#	   rrrrrrrrrrrrrrrrrr (goes into right contig)
	# This ended up with all kinds of base segments with negative
	# positions (since they were before the consensus.)

	if ( $right_cont_start < $right_cont_stop )
	{
		@right_bs = @{$right_contig->calculate_base_segments($right_cont_start, $right_cont_stop,$best_quality_reads)};
	}

	#There could be a situation in which we are tearing a contig like this:
	#	--------------------------------------------
	#	------------------------------
	#	-----------------------
	#				 ++++++++++++++++++++++++++
	#					++++++++++++++++++++++++
	#					   +++++++++++++++++++++
	# where the --- go into the new left contig and the
	# +++ go into the new right contig.
	# In this case, there is *no* region to the right of the torn
	# region.  The new right contig is formed entirely from the torn
	# region.  That is the reason for the if( ) below

	if ( $split_region_stop < $contig->length )
	{
		my $nBoundaryBaseSegment =
		  $contig->get_bs_index_at_position($split_region_stop +1, \@bs,  );

		my $rBoundaryBaseSegment = $bs[$nBoundaryBaseSegment];

		# create a new base segment which starts precisely at
		# stop + 1.  This may be the same as this
		# base segment, or it may be to the right of it.
		my $rNewBoundaryBaseSegment = {
			  type      => 'base_segment',
			  start => $split_region_stop + 1 + $offset,
			  stop => ($rBoundaryBaseSegment->{stop} + $offset),
			  read_name => $rBoundaryBaseSegment->{read_name},
		};

		push( @right_bs, $rNewBoundaryBaseSegment );

		for ( my $nBaseSeg = $nBoundaryBaseSegment + 1 ;
			  $nBaseSeg < @bs ;
			  $nBaseSeg++ )
		{
			my $rOldBaseSeg = $bs[$nBaseSeg];

			my $rNewBaseSeg = {
					type      => 'base_segment',
					start => $rOldBaseSeg->{start} + $offset,
					stop => ( $rOldBaseSeg->{stop} + $offset ),
					read_name => $rOldBaseSeg->{read_name},
			};

			push( @right_bs, $rNewBaseSeg );
		}
	}
	$right_contig->base_segments( \@right_bs);

	# do full checking once
	assert( $left_contig->is_bs_array_structure_ok( 1 ));
	assert( $right_contig->is_bs_array_structure_ok( 1 ));

	# modify consensus values in torn region in right contig

	my $clipped_start;
	my $clipped_stop;
	if (
		 calculate_intersection(
					 $right_cont_start,
					 $right_cont_stop,
					 1,
					 $right_contig->length,
					 \$clipped_start,
					 \$clipped_stop
		 )
	  )
	{
		Finishing::Assembly::Utility::recalculate_consensus_qualities_and_change($right_contig,
			$clipped_start,
			$clipped_stop,
			\@right_reads);
	}

	_transfer_consensus_tags($contig,$left_contig,$right_contig,$split_region_start,$split_region_stop, $offset);

	printf("tear complete.\n");    	

	return ($left_contig,$right_contig);		
}

sub _transfer_consensus_tags
{
	my ($contig, $left_contig, $right_contig, $start, $stop, $offset) = @_;
	my @left_tags,
	my @right_tags;
	foreach my $tag (@{$contig->tags})
	{

		if (    $tag->stop <= $start
			 && $tag->start < $stop )
		{

			# like this:
			# --------------------
			#           |              |
			# or this:
			# ---------
			#           |              |
			$tag->parent($tag->parent.'a');
			push( @left_tags, $tag );

		}
		elsif (    $start <= $tag->start
				&& $stop < $tag->stop )
		{

			# like this:
			#                     ----------------------------
			#              |            |
			# or this:
			#                                  --------------------
			#              |            |

			$tag->start($tag->start + $offset);
			$tag->stop($tag->stop + $offset);
			$tag->parent($tag->parent.'b');
			push( @right_tags, $tag );
		}
		else
		{

			# What is left is this:
			# ---------------------------------
			#              |            |
			# or this:
			#                ---------
			#              |            |

			# split the tag
			
            #my $right_tag = $tag->copy($tag);	No copy function found		
            my $right_tag = $tag;
			my $left_tag  = $tag;

			my $nOldConsPosEnd = $tag->stop;
			$left_tag->stop( min( $stop, $tag->stop ));
			$left_tag->parent($left_tag->parent .'a');

			$right_tag->start( max( $right_tag->start, $start ));
			$right_tag->start($right_tag->start + $offset);

			$right_tag->stop( $nOldConsPosEnd + $offset);

			$right_tag->parent($right_tag->parent . 'b');

			printf("had to split tag of type %s into two:  one from %d to %d in contig %s and the other from %d to %d in contig %s",
				$tag->type,                 $left_tag->start,
				$left_tag->stop,    $left_contig->name,
				$right_tag->start, $right_tag->stop,
				$right_contig->name
			);
            
            push @left_tags, $left_tag;
            push @right_tags, $right_tag;
		}
	}
	$left_contig->tags(\@left_tags);
	$right_contig->tags(\@right_tags);
}


#This split functions handles all fileio, PSE creation, etc.
sub read_data_and_split_contig
{
	my ($self, $data_source, $split_contig_name, $split_position, $no_gui, $out_file_name) = @_;
	
	my $cwd = getcwd;
	
	my $out_fh;
	my $bUsingDB=0;

	if(!($data_source =~ /\.ace/))
	{
		$bUsingDB=1
	}
	my $ace_object = $self->_get_ace_object($data_source, $split_contig_name);
	my $phd_object = $self->_get_phd_object($data_source);		
	
	my $rContig = $ace_object->get_contig($split_contig_name,1);

	my ($rLeftContig, $rRightContig) = $self->split($rContig, $phd_object, split_position => $split_position, no_gui => $no_gui);

	$ace_object->remove_contig($rContig->name);
	$ace_object->add_contig($rLeftContig);
	$ace_object->add_contig($rRightContig);

	if(!defined ($out_file_name))
	{	
		opendir THISDIR, ".";
		my @allFiles = readdir THISDIR;

		my $nMax = 1;
		for(my $i = 0;$i<@allFiles;$i++)
		{
			if($allFiles[$i] =~ /$data_source/)
			{
				my ( $Name, $Ext ) = $allFiles[$i] =~ /^(.+\.)(.+)$/;
				if(int($Ext)>$nMax)
				{
					$nMax=int($Ext);
				}
			}
		}
		$nMax++;
		$out_file_name = $data_source.".".$nMax;
	}
	my $file;
	if($bUsingDB)
	{
		$file = IO::String->new;
	}
	else
	{
		$file = IO::File->new(">$out_file_name");
	}
	$ace_object->write_file(output => $file);

#	if($bUsingDB)
#	{
#		$file->setpos(0);
#		my $split_contig;		
#
#		my $assembly_name = $data_source;
#		
#		my $assembly = GSC::Sequence::Assembly->get(
#   			sequence_item_name => $assembly_name,
#		);
#
#		$split_contig = GSC::Sequence::Contig->get(	sequence_item_name => "$assembly_name.$split_contig_name");		
#
#		my $ps = GSC::ProcessStep->get(process_to => 'split contigs'); # grab the 
#
#		my $pse = $ps->execute(
#    		contig => [$split_contig],
#			ace_fh => $file,
#    		assembly     => $assembly,    	
#			control_pse_id => 0,
#		);
#
#		unless ($pse){ # executes the pse, loads the new data into the 
#								# database
#    							# it didn't work for some reason;
#    		print $ps->error_message;
#		}
#
#		$split_contig->unlock;
#		
#
#		print "Commit of $assembly_name.$split_contig_name succeeded.\n\n";
#		App::DB->sync_database;
#		App::DB->commit;
#
#	}
}

sub complement
{	
	my ($self, $contig) = @_;
	#first, complement the consensus
	my $nConsLength = $contig->length;
	my $temp = reverse $contig->padded_base_string;
	$temp =~ tr/actgACTG/tgacTGAC/;
	$contig->padded_base_string($temp);
	$contig->complemented(!$contig->complemented);
	
	#next, reverse quality values for the consensus
	$contig->unpadded_base_quality([reverse @{$contig->unpadded_base_quality}]);

	#next complement reads and read positions for each read in the read array
	my $reads = $contig->assembled_reads;
    
	foreach my $read ($reads->all)
	{
        #my $temp = reverse $read->padded_base_string;
        my $temp = reverse $read->base_string;

		$temp =~ tr/actgACTG/tgacTGAC/;
        #$read->padded_base_string ($temp);
        $read->base_string($temp);

		$read->position(  $nConsLength - ($read->position + $read->length)+2);
		$read->complemented(!$read->complemented);
		
		my $ReadTags = $read->tags;
		next unless defined ($ReadTags);
		#complement read tags
		foreach my $tag (@{$ReadTags})
		{
			my $nReadLength = $read->length;
			#first change the offset, so that the start position and end
			#positions are now offsets from the end of the contig instead
			#of the beginning of the contig
			$tag->start( $nReadLength - $tag->start+1); 
			$tag->stop ($nReadLength - $tag->stop+1);
			#now swap start and end positions
			my $temp = $tag->stop;
			$tag->stop( $tag->start);
			$tag->start($temp); 	
		}
		$read->tags($ReadTags);
	}
    
	#now complement base segments
	my $base_segs = $contig->base_segments;
	foreach my $bs (@{$base_segs})
	{
		$bs->{start} = $nConsLength - $bs->{start}+1; 
		$bs->{stop} = $nConsLength - $bs->{stop}+1;
		my $temp = $bs->{stop};
		$bs->{stop} = $bs->{start};
		$bs->{start} = $temp;
	}

	@{$base_segs} = reverse @{$base_segs};
	$contig->base_segments($base_segs);
	
	

	#complement contig tags, if there are any for this particular contig
	my $ContigTags = $contig->tags;
	foreach my $tag (@{$ContigTags})
	{
		$tag->start( $nConsLength - $tag->start+1);
		$tag->stop ( $nConsLength - $tag->stop+1);

		#now swap start and stop positions
		my $temp = $tag->stop;
		$tag->stop ( $tag->start );
		$tag->start ( $temp);	
	}
	$contig->tags($ContigTags)if defined ($ContigTags);
    
    return $contig;
}

sub read_data_and_complement_contig
{
	my ($self, $data_source, $contig_name, %params) = @_;
	
	my $output_file_name = delete $params{output_file_name};
	my $using_db = 0;
	if(!($data_source =~ /\.ace/))
	{
		$using_db=1
	}

	my $ace_object = $self->_get_ace_object($data_source, $contig_name);
	my $phd_object = $self->_get_phd_object($data_source);
	
	my $contig = $ace_object->get_contig($contig_name);
	$self->complement($contig);
	$ace_object->add_contig($contig);

	my $file;
	if($using_db)
	{
		$file = IO::String->new;
	}
	else
	{
		$file = IO::File->new(">$output_file_name");
	}

	$ace_object->write_file(output => $file);
	
#	if($using_db)
#	{
#		$file->setpos(0);
#		my $merge_contig;		
#
#		my $assembly_name = $data_source;
#		
#		my $assembly = GSC::Sequence::Assembly->get(
#   			sequence_item_name => $assembly_name,
#		);		
#		
#		my $db_contig = GSC::Sequence::Contig->get(	sequence_item_name => "$assembly_name.$contig_name");		
#		my $ps = GSC::ProcessStep->get(process_to => 'complement contig'); # grab the process step
#
#		my $pse = $ps->execute(
#    		old_contig => [$db_contig],
#			ace_fh => $file,
#    		assembly => $assembly,    	
#			control_pse_id => 0,
#		);
#
#		unless ($pse){ # executes the pse, loads the new data into the 
#								# database
#    							# it didn't work for some reason;
#    		print $ps->error_message;
#		}
#
#		$db_contig->unlock;
#
#		print "Commit of $assembly_name.$contig_name succeeded.\n\n";
#		App::DB->sync_database;
#		App::DB->commit;
#
#	}

}

sub read_data_and_merge_contigs
{
	my ($self, $data_source, $left_contig_name, $right_contig_name, %params) = @_;
	
	my $bUsingDB=0;

	if(!($data_source =~ /\.ace/))
	{
		$bUsingDB=1
	}
	
	my $output_file_name = delete $params{output_file_name};
	my $ace_object = $self->_get_ace_object($data_source, $left_contig_name, $right_contig_name);
	my $phd_object = $self->_get_phd_object($data_source);

	my $left_contig = $ace_object->get_contig($left_contig_name);
	my $right_contig = $ace_object->get_contig($right_contig_name);
	my $merge_contig = $self->merge($left_contig, $right_contig, $phd_object, %params);

	$ace_object->remove_contig($left_contig_name);
	$ace_object->remove_contig($right_contig_name);
	$ace_object->add_contig($merge_contig);
	
	my $file;
	if($bUsingDB)
	{
		$file = IO::String->new;
	}
	else
	{
		$file = IO::File->new(">$output_file_name");
	}

	$ace_object->write_file(output => $file);
	
#	if($bUsingDB)
#	{
#		$file->setpos(0);
#		my $merge_contig;		
#
#		my $assembly_name = $data_source;
#		
#		my $assembly = GSC::Sequence::Assembly->get(
#   			sequence_item_name => $assembly_name,
#		);
#		
#		my $db_left_contig = GSC::Sequence::Contig->get(	sequence_item_name => "$assembly_name.$left_contig_name");		
#		my $db_right_contig = GSC::Sequence::Contig->get(	sequence_item_name => "$assembly_name.$right_contig_name");		
#		my $ps = GSC::ProcessStep->get(process_to => 'merge contigs'); # grab the 
#
#		my $pse = $ps->execute(
#    		old_contigs => [$left_contig, $right_contig],
#			ace_fh => $file,
#    		assembly => $assembly,    	
#			control_pse_id => 0,
#		);
#
#		unless ($pse){ # executes the pse, loads the new data into the 
#								# database
#    							# it didn't work for some reason;
#    		print $ps->error_message;
#		}
#
#		$db_left_contig->unlock;
#		$db_right_contig->unlock;
#
#		print "Commit of $assembly_name.$left_contig_name and $assembly_name.$right_contig_name succeeded.\n\n";
#		App::DB->sync_database;
#		App::DB->commit;
#
#	}


}

sub remove_reads
{
	my ($self, $ori_contig, $phd_object, %params) = @_;		

	my %reads = map {$_->name, $_} $ori_contig->assembled_reads->all;
	my $remove_reads = $params{remove_reads};
    
	foreach my $read_name (@{$remove_reads}) 
	{ 
		if(delete $reads{$read_name})
		{
			print "remove $read_name from ".$ori_contig->name."\n";
		}
		else
		{
			warn "Could not find $read_name in ".$ori_contig->name."\n";
		}
	}
	if((scalar keys %reads) == 0)
	{
		print $ori_contig->name." has no reads left.\n";
		return undef;
	}	
    #my $region = $ori_contig->get_contiguous_reads([values %reads],$ori_contig->tags);
    #if(!defined $region) {print $ori_contig->name."\n"; return undef;}

    my $src_fac = Finishing::Assembly::Factory->connect('source'); 
    my $contig  = $src_fac->copy_contig($ori_contig);

    #$contig->assembled_reads($region->{reads});
    $contig->assembled_reads([values %reads]);
    #%reads = map{$_->name, $_}@{$region->{reads}};
    
    #$contig->tags($region->{tags});
    $contig->tags($ori_contig->tags);
	# at this point, all reads in the contig have been put into 2 lists:
	# a list of reads that the user wants to go into the 'left' contig,
	# and a list of reads that the user wants to go into the 'right contig.
	# (Note that these contigs themselves may fall apart into more contigs.)

	# Find the torn region--To the left and right of this region, any column
	# will be precisely the same as it was before the tear:  the same reads
	#  be in that column.  Within the tear, the reads will be different,
	# and thus we will have to recompute the consensus, Golden path, and
	# base segment.

	my $start = 1;
    #my $stop = $region->{length};
    my $stop  = $ori_contig->length;

=cut

	my @read_names = ();
    
	foreach my $read_name (keys %reads)
	{
		my ($version) = ($reads{$read_name}->phd_file =~ /\.(\d+)$/); 
		$read_name .= "-$version";
        push @read_names, $read_name;
	}	
	my @reads = GSC::Sequence::Item->get(sequence_item_name => \@read_names);
	foreach my $read ( @reads)
	{
        my $qual = $read->sequence_quality_arrayref;
        unless ($qual) {
            print $read->sequence_item_name."has no sequence_quality_arrayref\n";
            next;
        }
		my @quality = @{$read->sequence_quality_arrayref};
		die "Couldln't load $read->sequence_item_name quality.\n" if ((scalar @quality) == 0);
		my ($read_name) = $read->sequence_item_name =~ /(.+)-\d+/;
		die "Couldn't resolve read name\n" if ! exists $reads{$read_name}; 
		@quality = reverse @quality if ($reads{$read_name}->complemented);#print $read_name."\n";		
		$reads{$read_name}->unpadded_base_quality( \@quality);				
	}

=cut

	foreach( values %reads)
	{
		my $phd = $phd_object->get_phd($_->phd_file);
        unless ($phd) {
            print 'phd_file Not Found for '.$_->name."\n";
            exit;
        }
		my @quality = @{$phd->unpadded_base_quality};
		@quality = reverse @quality if ($_->complemented);
		$_->unpadded_base_quality( \@quality);				
	}
	print "Done loading read quality for ".$contig->name."\n";

	exit if ( !check_contiguous( [values %reads], $start, $stop ) );	
	
	# point of no return
	# At this point, there is no longer any error checking.
	# Errors indicate a program bug.

	# in the torn region, figure out the consensus base and quality
	# The fast method of doing this is to go through each read once,
	# and mark two arrays--the highest quality at a cons pos and the
	# read that attained that highest quality.  I will ignore the
	# "aligned" clipping at this point, since a read could be high
	# quality unaligned precisely because there is a mis-join, which
	# is the reason the user is tearing the contig.  Thus the correct
	# consensus could be this unaligned high quality region.
	
	
	assert( $start <= $stop );
	my ($best_quality_reads, $best_quality) = $contig->calculate_consensus($start,$stop,[values %reads]);
		
	# transferring base segments

	# which is the last base segment that can be transferred?

	# There could be a situation in which we are tearing a contig like this:
	#			cccccccccccccccccccccccccc
	#	---------------------------
	#		 +++++++++++++++++++++++++++++++++++
	#
	# where the --- go into the new left contig and the
	# +++ go into the new right contig.
	# In this case, the torn region starts before the consensus
	# and thus there is no base segments before the torn region.
	# Thus no transferring of base segments before the torn region.

	$contig->base_segments($contig->calculate_base_segments($start, $stop,$best_quality_reads));
		
	# now recalculate the consensus quality values for the new left contig
	# (The reason for the intersect is that there may be reads that
	# stick out on the left to negative consensus positions and some
	# of these may be in the right contig and thus the start_
	# may be to the left of where the new left contig starts.)

	Finishing::Assembly::Utility::recalculate_consensus_qualities_and_change($contig, $start, $stop, [values %reads]);
	
	# That completes the base segments for the new left contig
	# do full checking once
	assert( $contig->is_bs_array_structure_ok( 1 ));

	printf("read removal complete.\n");    	

	return $contig;		
}

sub remove_read
{
    my ($self, $ori_contig, $phd_object, %params) = @_;		

    my %reads = map{$_->name, $_}$ori_contig->assembled_reads->all;
    my $read = delete $reads{$params{remove_read}};

    my $src_fac = Finishing::Assembly::Factory->connect('source'); 
    my $contig  = $src_fac->copy_contig($ori_contig);
    $contig->assembled_reads([values %reads]);
    $contig->tags($ori_contig->tags);
    
	# at this point, all reads in the contig have been put into 2 lists:
	# a list of reads that the user wants to go into the 'left' contig,
	# and a list of reads that the user wants to go into the 'right contig.
	# (Note that these contigs themselves may fall apart into more contigs.)

	# Find the torn region--To the left and right of this region, any column
	# will be precisely the same as it was before the tear:  the same reads
	#  be in that column.  Within the tear, the reads will be different,
	# and thus we will have to recompute the consensus, Golden path, and
	# base segment.

	my $start = $read->align_start;
	my $stop = $read->align_stop;	
	# now lets find if things fall into more than 2 contigs in the torn region
	# To do this, let's make two little lists:  The left reads in the
	# torn region, and the right reads in the torn region.
	my $overlap_reads = $contig->get_overlapping_reads($start, $stop, [values %reads]);

	foreach (@$overlap_reads)
	{
		my $phd = $phd_object->get_phd($_->phd_file);
		my @quality = @{$phd->unpadded_base_quality};
		@quality = reverse @quality if ($_->complemented);
		$_->unpadded_base_quality( \@quality);				
	}
	# Now let's check (for the left contig and the right contig)
	# that there is a base for every position
	# (unaligned or aligned) for this entire region
	my $cont_start = max( $start,  1 );
	my $cont_stop   = min( $stop, $contig->length );
	exit if ( !check_contiguous( $overlap_reads, $cont_start, $cont_stop ) );	
	
	# point of no return
	# At this point, there is no longer any error checking.
	# Errors indicate a program bug.

	# in the torn region, figure out the consensus base and quality
	# The fast method of doing this is to go through each read once,
	# and mark two arrays--the highest quality at a cons pos and the
	# read that attained that highest quality.  I will ignore the
	# "aligned" clipping at this point, since a read could be high
	# quality unaligned precisely because there is a mis-join, which
	# is the reason the user is tearing the contig.  Thus the correct
	# consensus could be this unaligned high quality region.
	
	
	assert( $cont_start <= $cont_stop );
	my ($best_quality_reads, $best_quality) = $contig->recalculate_consensus($cont_start,$cont_stop,$overlap_reads);
		
	# transferring base segments

	# which is the last base segment that can be transferred?

	# There could be a situation in which we are tearing a contig like this:
	#			cccccccccccccccccccccccccc
	#	---------------------------
	#		 +++++++++++++++++++++++++++++++++++
	#
	# where the --- go into the new left contig and the
	# +++ go into the new right contig.
	# In this case, the torn region starts before the consensus
	# and thus there is no base segments before the torn region.
	# Thus no transferring of base segments before the torn region.

	my @bs = @{$ori_contig->base_segments};
	my @bs_new = @{$ori_contig->calculate_base_segments($cont_start, $cont_stop,$best_quality_reads)};
	my $bs_start = $ori_contig->get_bs_index_at_position($cont_start);
	my $bs_end = $ori_contig->get_bs_index_at_position($cont_stop);
	$bs[$bs_start-1]->{stop} = $cont_start-1 if($cont_start>1);
	$bs[$bs_end+1]->{start} = $cont_stop+1 if($cont_stop<$ori_contig->length);
	splice(@bs, $bs_start,$bs_end-$bs_start+1,@bs_new);
	$contig->base_segments( \@bs);
	
	# now recalculate the consensus quality values for the new left contig
	# (The reason for the intersect is that there may be reads that
	# stick out on the left to negative consensus positions and some
	# of these may be in the right contig and thus the start_
	# may be to the left of where the new left contig starts.)

	my $nConsPosLeftToRecalculate;
	my $nConsPosRightToRecalculate;

	if (
	    calculate_intersection(
				   1,                      $contig->length,
				   $start,             $stop,
				   \$nConsPosLeftToRecalculate, \$nConsPosRightToRecalculate
				   )
	    )
	{
	    Finishing::Assembly::Utility::recalculate_consensus_qualities_and_change($contig, $nConsPosLeftToRecalculate, $nConsPosRightToRecalculate, $overlap_reads);
	}

	# That completes the base segments for the new left contig


	# do full checking once
	assert( $contig->is_bs_array_structure_ok( 1 ));

	printf("read removal complete.\n");    	

	return $contig;		
}

sub replace_xns
{
    my ($self, $contig) = @_;

    #GET A LIST OF N/X POSITIONS IN CONTIG

    my $xn_positions = $contig->padded_base_string_xn_positions;

    print "No x/ns found in contig: " . $contig->name . "\n" and return $contig
	unless scalar @{$xn_positions} > 0;

    #GET THE ACTUAL BASE AT THAT POSITION FROM EACH READ
    #ALIGNED UNDER THE N/X POSITION

    my $reads = $contig->assembled_reads;
    my $base_counts = {};
    my $replace_xns = {};

    foreach my $pos (@$xn_positions)
    {
	foreach my $read ($reads->all)
	{
	    my @read_positions = ($read->start .. $read->stop);
	    next unless grep (/^$pos$/, @read_positions);
	    my $read_position = $pos - $read->start;
	    my $base = $read->get_base_at_position ($read_position);

	    #CREATE HASH OF BASE AND NUMBER OF TIMES THAT BASE IS
	    #PRESENT 

	    next if $base =~ /[nx]/i;
	    $base_counts->{$base}++;
	}

	my @base_candidates;
	my $prev_count = 0;
	my $curr_count = 0;

	#GRAB THE HIGHEST REPRESENTED BASE OR BASES
	#IF THERE'S A TIE ONE WILL BE CHOSE RANDOMLY

	foreach my $base (sort {$base_counts->{$b} <=> $base_counts->{$a}} keys %$base_counts)
	{
	    $curr_count = $base_counts->{$base};
	    next if $curr_count < $prev_count;
	    push @base_candidates, $base if $curr_count >= $prev_count;
	    $prev_count = $curr_count;
	}

	my $new_base;

	if (scalar @base_candidates == 0)
	{
	    #NO BASE CALLED BY ANY READ SO PICK ONE RANDOMLY
	    #OR ALL READS CALL AN X OR N
	    @base_candidates = qw/ a c g t /;
	    $new_base = @base_candidates[int rand (scalar @base_candidates)];
	}
	elsif (scalar @base_candidates == 1)
	{
	    #THERE IS ONE BASE THAT HAS THE HIGHEST REPRESENTATION
	    $new_base = @base_candidates[0];
	}
	else
	{
	    #THERE ARE MULTIPLE READS WITH SAME HIGHEST REPRESENTATION
	    #PICK ONE RANDOMLY
	    $new_base = @base_candidates[int rand (scalar @base_candidates)];
	}

	#CREATE A HASH OF POSITION AND NEW BASE

	$replace_xns->{$pos} = $new_base;
    }

    #REPLACE THE BASE
    
    my $base_position = 0;

    my $new_base_string;
    my @new_quals;
    
    foreach my $base (split ('', $contig->padded_base_string))
    {
	$base_position++;

	#MAKE SURE BASE AND QUALITY ALIGN CORRECTLY

	if ($base =~ /^[xn]$/i)
	{
	    print "Unable to get replacement base for x/n at position $base_position\n"
		and return unless exists $replace_xns->{$base_position};

	    $new_base_string .= $replace_xns->{$base_position};

	    push @new_quals, 0;

	    next;
	}

	$new_base_string .= $base;

	push @new_quals, ${$contig->padded_base_quality}[$base_position - 1];

	#SET REPLACED QUAL TO ZERO

    }

    #REPLACE BASE STRING AND QUALITIES 
    
    $contig->padded_base_string ($new_base_string);
    $contig->padded_base_quality (\@new_quals);

    return $contig;
}


=pod

=head1 NAME

ContigTools - Object oriented contig toolkit

=head1 SYNOPSIS

 my ($left_contig, $right_contig) = ContigTools->split($contig, $phd_obj, split_position => $position);

 my $merge_contig = ContigTools->merge($left_contig, $right_contig, $phd_obj);
    
=head1 DESCRIPTION

Finishing::Assembly::ContigTools performs operation on objects of type contig.  So far the functions supported are merging and splitting.

=head1 METHODS

=head1  merge 

my $merge_contig = ContigTools->merge($left_contig, $right_contig, $phd_obj);

merges the ends of the two contigs together, performing an alignment in the process.

=head1 split 

my ($left_contig, $right_contig) = ContigTools->split($contig, $phd_obj, split_position => $position);
    
splits the contig give at the specified padded position coordinates.

=head1 AUTHOR

This module was written and is maintained by Jon Schindler <jschindl@wugsc.wustl.edu>.
    
=cut

1;
#$HeadURL: svn+ssh://svn/srv/svn/gscpan/perl_modules/branches/jschindl/Assembly/ContigTools.pm $
#$Id: ContigTools.pm 22733 2007-03-01 02:25:45Z jschindl $
