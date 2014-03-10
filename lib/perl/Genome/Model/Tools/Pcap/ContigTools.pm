package Genome::Model::Tools::Pcap::ContigTools;

our $VERSION = 0.01;

use strict;

use Genome::Model::Tools::Pcap::Transform;
use Carp::Assert;
use Cwd;
use Bio::Seq;
use Bio::Tools::dpAlign;
use Genome::Model::Tools::Pcap::Utility;
use Genome::Model::Tools::Pcap::Ace::Writer;
use Genome::Model::Tools::Pcap::Ace::Reader;
use Genome::Model::Tools::Pcap::Ace;
use List::Util qw(min max);

our $initialized = 0;
sub init_gtk {
    return if $initialized;
    eval {
        require Gtk;
        require Gtk::ButtonCrate::Radio;        
    };
    die $@ if $@;
    $initialized = 1;
}

sub new 
{
    croak("no class given, quitting") if @_ < 1;
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

    #lets grab the last 900 bases of the left contig, and the first 900 bases
    # of the right contig.  We only need to do this if they are longer than
    # 900 bases, which is why I'm checking

    if(length($left_string) >= 900)
    {
        $left_string = substr($left_string, length($left_string)-900, 900);	
    }

    if(length($right_string) >= 900)
    {
        $right_string = substr($right_string, 0, 900);	
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
    $left_contig->{align} = {};
    $left_contig->{align}{start}= $out->{_start_end_lists}{$left_contig->name}[0]{start};
    $left_contig->{align}{end} = $out->{_start_end_lists}{$left_contig->name}[0]{end};

    #appears in the left contig to the one produced by bioperl
    #we use it later when adding these pads to the reads	
    $left_contig->{align}{local_align_seq} = $left_contig->{align}{align_seq} = $out->{_start_end_lists}{$left_contig->name}[0]{seq};
    $right_contig->{align}{local_align_seq} = $right_contig->{align}{align_seq} = $out->{_start_end_lists}{$right_contig->name}[0]{seq};
    $left_contig->{align}{align_seq} =~ tr/rR/**/;
    $right_contig->{align}{align_seq} =~ tr/rR/**/;
    if( $left_contig->length > 900)
    {
        $left_contig->{align}{start} = $left_contig->length-900+$left_contig->{align}{start};
        $left_contig->{align}{end} = $left_contig->length-900+$left_contig->{align}{end};
    }
    else
    {
        $left_contig->{align}{start} =	$left_contig->{align}{start};
        $left_contig->{align}{end} = $left_contig->{align}{end};
    }
    $right_contig->{align}{start}= $out->{_start_end_lists}{$right_contig->name}[0]{start};
    $right_contig->{align}{end} = $out->{_start_end_lists}{$right_contig->name}[0]{end};
    my $left_pre = substr( $left_contig->padded_base_string, 0, $left_contig->{align}{start}-1 );
    my $left_post = substr( $left_contig->padded_base_string, $left_contig->{align}{end}, $left_contig->length-$left_contig->{align}{end});
    $left_contig->{align}{merge_seq} = $left_pre.$left_contig->{align}{align_seq};
    $left_contig->{align}{align_seq} = $left_pre.$left_contig->{align}{align_seq}.$left_post;	
    my $right_pre = substr( $right_contig->padded_base_string, 0, $right_contig->{align}{start}-1 );
    my $right_post = substr( $right_contig->padded_base_string, $right_contig->{align}{end}, $right_contig->length-$right_contig->{align}{end});
    $right_contig->{align}{align_seq} = $right_pre.$right_contig->{align}{align_seq}.$right_post;
    $left_contig->{align}{merge_seq} .= $right_post;
    $left_contig->{align}{merge_seq} =~ tr/-/*/;


    my $left_transform = Genome::Model::Tools::Pcap::Transform->new($left_contig->{align}{align_seq}, '-');#this creates a transform for the merge region as it originally 
    my $right_transform = Genome::Model::Tools::Pcap::Transform->new($right_contig->{align}{align_seq}, '-');#this creates a transform for the merge region as it originally 
    $right_transform->_offset( $left_contig->{align}{start}-$right_contig->{align}{start});
    #may want to check if offset is valid
    $left_contig->{align}{transform} = $left_transform;
    $right_contig->{align}{transform} = $right_transform;
}

sub _clip_read_after_position
{
    my ($read, $transform,$position) = @_;
    my $clip_position = $read->get_child_position_from_parent_position($position);
    $read->align_clip_start(min(max($read->align_clip_start,$clip_position),$read->align_end));
    $read->qual_clip_start(min(max($read->qual_clip_start,$clip_position),$read->align_end));
    $read->qual_clip_end(min(max($read->qual_clip_end,$clip_position),$read->align_end));
    $read->align_clip_end(min(max($read->align_clip_end,$clip_position),$read->align_end));

    $read->align_clip_start($transform->get_pad_position($read->align_clip_start));		
    $read->align_clip_end($transform->get_pad_position($read->align_clip_end));		
    $read->qual_clip_start($transform->get_pad_position($read->qual_clip_start));		
    $read->qual_clip_end($transform->get_pad_position($read->qual_clip_end));
    if($read->align_clip_start>=$read->align_clip_end||
        $read->qual_clip_start>=$read->qual_clip_end )
    {
        $read->align_clip_start($read->align_clip_end($read->qual_clip_start($read->qual_clip_end(1))));	
    }
}

sub _clip_read_before_position
{
    my ($read, $transform,$position) = @_;
    my $clip_position = $read->get_child_position_from_parent_position($position);
    $read->align_clip_start(max(min($read->align_clip_start,$clip_position)),1);
    $read->qual_clip_start(max(min($read->qual_clip_start,$clip_position)),1);
    $read->qual_clip_end(max(min($read->qual_clip_end,$clip_position)),1);
    $read->align_clip_end(max(min($read->align_clip_end,$clip_position)),1);

    $read->align_clip_start($transform->get_pad_position($read->align_clip_start));		
    $read->align_clip_end($transform->get_pad_position($read->align_clip_end));		
    $read->qual_clip_start($transform->get_pad_position($read->qual_clip_start));		
    $read->qual_clip_end($transform->get_pad_position($read->qual_clip_end));
    if($read->align_clip_start>=$read->align_clip_end||
        $read->qual_clip_start>=$read->qual_clip_end)
    {
        $read->align_clip_start($read->align_clip_end($read->qual_clip_start($read->qual_clip_end(1))));	
    }

}

sub _create_merge_tag
{
    my ($contig_name,$left_end_of_merge, $right_end_of_merge) = @_;
    my @temptime = localtime;
    my $year = sprintf("%02d", $temptime[5] % 100);
    my $tempstring = sprintf("%02d%02d%02d:%02d%02d%02d", $year, ($temptime[4]+1), $temptime[3], 
        $temptime[2], $temptime[1], $temptime[0] );

    return Genome::Model::Tools::Pcap::Tag->new(type => 'comment', , parent =>
        $contig_name, start => $left_end_of_merge, stop => $right_end_of_merge, date => $tempstring,
        source => "cmt", no_trans => "NoTrans", scope => "ACE");
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
        $tag->start ( $right_contig->{align}{transform}->get_pad_position($tag->start, 0));
        $tag->stop ( $right_contig->{align}{transform}->get_pad_position($tag->stop, 0));							
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
        $stats{start} = $left_contig->{align}{transform}->get_pad_position($left_contig->{align}{start});
        $stats{stop} = $left_contig->{align}{transform}->get_pad_position($left_contig->{align}{end});		
        $stats{left_contig} = $left_contig;
        $stats{right_contig} = $right_contig;
    }

    $stats{mismatch} =0;
    $stats{length} = length($left_contig->{align}{local_align_seq});
    for(my $pos=0;$pos<$stats{length};$pos++)
    {
        if(lc(substr($left_contig->{align}{local_align_seq},$pos,1)) ne lc(substr($right_contig->{align}{local_align_seq}, $pos, 1)))
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
    for(my $pos = $left_contig->{align}{start};$pos <= $left_contig->{align}{end};$pos++)
    {
        push @qarray1, $left_contig->padded_base_quality->[$pos-1];
    }
    my @qarray2;
    for(my $pos = $right_contig->{align}{start};$pos <= $right_contig->{align}{end};$pos++)
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
            if(lc(substr($left_contig->{align}{local_align_seq},$pos,1)) ne lc(substr($right_contig->{align}{local_align_seq}, $pos, 1)))
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
            if(lc(substr($left_contig->{align}{local_align_seq},$pos,1)) ne lc(substr($right_contig->{align}{local_align_seq}, $pos, 1)))
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

    print $fh "Merge sequence is:\nLeft Contig:  $left_contig->{align}{local_align_seq}\n";
    print $fh "Right Contig: $right_contig->{align}{local_align_seq}\n";

    print $fh "Discrepancy:  $stats->{match_string}\n\n";
    print $fh "Statistics Info:\n";
    print $fh "Length of merge region:               $stats->{length}\n";
    print $fh "Number of discrepancies:              $stats->{mismatch}\n";
    printf $fh ("Percent Identity:                     %2.2f\%\n", $stats->{percent_identity});
    print $fh "Number of high quality bases:         $stats->{hqcount}\n";
    print $fh "Number of high quality discrepancies: $stats->{hqmismatch}\n";

    printf $fh ("High Quality Base Percent Identity:   %2.2f\%\n", $stats->{hq_percent_identity});

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
        my $temp1 = $left_contig->{align}{transform}->get_pad_position($left_contig->{align}{start});
        my $temp2 = $left_contig->{align}{transform}->get_pad_position($left_contig->{align}{end});
        print "Merging ",$left_contig->name, " and ", $right_contig->name, " from consensus position $temp1 to $temp2. \n";
    }
    print "Merge sequence is:\nLeft Contig:  $left_contig->{align}{local_align_seq}\n";
    print "Right Contig: $right_contig->{align}{local_align_seq}\n";
    my $mismatch =0;
    my $length = length($left_contig->{align}{local_align_seq});
    for(my $pos=0;$pos<$length;$pos++)
    {
        if(lc(substr($left_contig->{align}{local_align_seq},$pos,1)) ne lc(substr($right_contig->{align}{local_align_seq}, $pos, 1)))
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
    for(my $pos = $left_contig->{align}{start};$pos <= $left_contig->{align}{end};$pos++)
    {
        push @qarray1, $left_contig->padded_base_quality->[$pos-1];
    }
    my @qarray2;
    for(my $pos = $right_contig->{align}{start};$pos <= $right_contig->{align}{end};$pos++)
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
            if(lc(substr($left_contig->{align}{local_align_seq},$pos,1)) ne lc(substr($right_contig->{align}{local_align_seq}, $pos, 1)))
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
            if(lc(substr($left_contig->{align}{local_align_seq},$pos,1)) ne lc(substr($right_contig->{align}{local_align_seq}, $pos, 1)))
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
    my $templength = length( $left_contig->{align}{local_align_seq} );
    print "Length of merge region:               $templength\n";
    print "Number of discrepancies:              $mismatch\n";
    printf ("Percent Identity:                     %2.2f\%\n", $percent_identity);
    print "Number of high quality bases:         $hqcount\n";
    print "Number of high quality discrepancies: $hqmismatch\n";

    printf ("High Quality Base Percent Identity:   %2.2f\%\n", $nHqPercentIdentity);
}

sub merge
{
    my ($self, $left_contig, $right_contig, $phd_object, %params) = @_;

    my $nGapOpenPenalty=$params{gap_open_penalty} || 1;
    my $nGapExtPenalty=$params{gap_ext_penalty} || 1;
    my $nGGapExtPenalty=$params{glob_gap_open_penalty} || 15;
    my $nGGapOpenPenalty=$params{glob_gap_open_penalty} || 15;
    my $bUseGlobalAlign=$params{use_global_align} || 0;
    my $bQuiet=$params{quiet} || 0;	

    my %left_reads = %{$left_contig->reads};
    my %right_reads = %{$right_contig->reads};
    my %left_overlap_reads;
    my %right_overlap_reads;
    my @base_segs;		
    assert($left_contig->is_bs_array_structure_ok(1));
    assert($right_contig->is_bs_array_structure_ok(1));

    _align_contigs($left_contig,$right_contig,1,1);		
    #print_stats($left_contig, $right_contig);
    my $stats = compute_stats($left_contig, $right_contig);
    print_stats2($stats, $params{statsfh});
    if(!check_stats($stats,$params{cutoffs})){ return };	

    #take the left contig up to the last 900 bases and repad it.

    #here is where we merge the consensus for the two contigs, basically, we know that the locally aligned region
    #will be the same for the left and right contigs.  So, we just add the consensus in the right contig that occurs
    #after the locally aligned region to the left contig.

    my $merge_contig = Genome::Model::Tools::Pcap::Contig->new;
    $merge_contig->name( $left_contig->name );
    $merge_contig->{align} = {};
    $merge_contig->{align}{start} = $left_contig->{align}{transform}->get_pad_position($left_contig->{align}{start});
    $merge_contig->{align}{end} = $left_contig->{align}{transform}->get_pad_position($left_contig->{align}{end});
    $merge_contig->complemented( $left_contig->complemented );
    $merge_contig->padded_base_string($left_contig->{align}{merge_seq});	

    foreach(values %left_reads)
    {

        if($_->align_end > $left_contig->{align}{start} )
        {
            my $new_seq = $left_contig->{align}{transform}->pad_string_partial($_->padded_base_string, $_->position-1, '-');
            my $transform = Genome::Model::Tools::Pcap::Transform->new($new_seq, '-');
            $new_seq =~ tr/-/*/;
            $_->padded_base_string( $new_seq);


            # set clipping to before the end of the merge region if necessary
            _clip_read_before_position($_,$transform,$left_contig->{align}{end});
            $_->{padded_base_count} = length($_->padded_base_string);
            $_->{fromcontig} = "left";
            $left_overlap_reads{$_->name()} = $_;
        }
        $_->position ($left_contig->{align}{transform}->get_pad_position($_->position));
    }

    foreach(values %right_reads)
    {
        if($_->position <= $right_contig->{align}{end} )
        {
            my $tempoffset = 0;#$right_contig->{align}{transform}->_offset;
            my $new_seq = $right_contig->{align}{transform}->pad_string_partial($_->padded_base_string, $_->position-$tempoffset-1, '-');
            my $transform = Genome::Model::Tools::Pcap::Transform->new($new_seq, '-');
            $new_seq =~ tr/-/*/;
            $_->padded_base_string( $new_seq);


            # set clipping to before the end of the merge region if necessary
            _clip_read_after_position($_,$transform,$right_contig->{align}{start});

            $_->{padded_base_count} = length($_->padded_base_string);
            $_->{fromcontig} = "right";
            $right_overlap_reads{$_->name()} = $_;
        }
        $_->position ($right_contig->{align}{transform}->get_pad_position($_->position));
    }

    foreach(values %left_overlap_reads)
    {
        # first, take the read up to the locally aligned region, and copy it
        # exactly as it is


    }

    foreach(values %right_overlap_reads)
    {
        # first, take the read up to the locally aligned region, and copy it
        # exactly as it is

    }

    #go ahead and add the reads together
    #add the left reads to the list of merge reads
    my %overlap_reads;
    foreach( values %left_overlap_reads,values %right_overlap_reads)
    {

        my $phd = $phd_object->get_phd($_->phd_file);
        my @quality = @{$phd->unpadded_base_quality};
        @quality = reverse @quality if ($_->complemented);
        $_->unpadded_base_quality( \@quality);
        $overlap_reads{$_->name} = $_;		
    }
    $merge_contig->reads({%left_reads,%right_reads});

    #now that we have our list of reads in the merge region, lets open the quality values
    #for these reads so that we can recalculate the quality values for the consensus

    my $left_end_of_merge = $merge_contig->{align}{start};
    my $right_end_of_merge = $merge_contig->{align}{end};

    # reset the consensus base and quality arrays to refigure them
    # base on the combined reads from both the left and right contig.  Note that it is
    # possible (likely) that the consensus bases and quality in the
    # merged region will be different for the new contig.  

    my ($best_quality_reads,$best_quality) = $merge_contig->recalculate_consensus($merge_contig->{align}{start},
        $merge_contig->{align}{end},
        \%overlap_reads);

    # append base segments for the left contig before the merge region
    my $bs_position = $left_contig->get_bs_index_at_position( $left_contig->{align}{start});
    push( @base_segs, @{$left_contig->base_segments} );
    #don't add BS where there aren't any, make sure the array is at least one unit long
    if($base_segs[$bs_position]{start_pos}<$merge_contig->{align}{start})
    {
        #in this case, the base segment at this position starts before the merge region
        #clip the end pos to the base before the merge region
        $base_segs[$bs_position]{end_pos} = $merge_contig->{align}{start}-1;	
        $#base_segs = $bs_position;
    }
    else
    {
        $#base_segs = $bs_position-1;
    }
    #assert( $merge_contig->is_bs_array_ok(0) );			

    my @bs_temp = @{$merge_contig->calculate_base_segments($merge_contig->{align}{start},$merge_contig->{align}{end},$best_quality_reads)};
    push(@base_segs, @bs_temp);

    #assert( $merge_contig->get_bs_array_structure_ok(0) );			
    # get the Base Segments that start after the merge region, and recalculate their start and end positions,
    # then push them onto the base seg array (after the merge region)
    $bs_position = $right_contig->get_bs_index_at_position( $right_contig->{align}{end}+1);
    @bs_temp = @{$right_contig->base_segments};
    @bs_temp = splice ( @bs_temp, $bs_position, @bs_temp - $bs_position);
    if(@bs_temp>0&&$bs_temp[0]->{end_pos}>=($right_contig->{align}{end}+1))
    {
        foreach my $bs (@bs_temp)
        {
            $bs->{start_pos} = $right_contig->{align}{transform}->get_pad_position($bs->{start_pos});
            $bs->{end_pos} = $right_contig->{align}{transform}->get_pad_position($bs->{end_pos});
        }
        #clip the first BS of aBSTemp so that it starts AFTER the merge
        #region!
        if($bs_temp[0]{start_pos}<=$base_segs[@base_segs-1]{end_pos})
        {
            $bs_temp[0]{start_pos} = $base_segs[@base_segs-1]{end_pos}+1;
        }

        push(@base_segs, @bs_temp);
        #$merge_contig->base_segments(\@base_segs);
        #assert( $merge_contig->is_bs_array_structure_ok(1)) ;

    }	
    $merge_contig->base_segments(\@base_segs);
    assert( $merge_contig->is_bs_array_structure_ok(1)) ;

    # next, splice in the newly computed quality values above to the left contig, then we will append the quality

    my @cons_quality = @{$left_contig->padded_base_quality};
    $#cons_quality = $left_contig->{align}{start}-2;

    my @temp = @{$best_quality};
    splice(@temp,0,$merge_contig->{align}{start});
    push (@cons_quality, @temp);

    my @right_cons_quality = @{$right_contig->padded_base_quality};	
    splice(@right_cons_quality,0,$right_contig->{align}{end});
    push (@cons_quality,@right_cons_quality);
    $merge_contig->padded_base_quality(\@cons_quality);
    # now recalculate consensus qualities

    Genome::Model::Tools::Pcap::Utility::recalculate_consensus_qualities_and_change( $merge_contig, $merge_contig->{align}{start}, $merge_contig->{align}{end}, [values %overlap_reads]  );

    #transfer left contig tags and translated right contig tags
    my @merge_contig_tags = (@{$left_contig->tags}, _translate_tag_positions($right_contig,$merge_contig));
    push @merge_contig_tags, _create_merge_tag($merge_contig->name,$merge_contig->{align}{start},$merge_contig->{align}{end});	
    _set_tags_parent_name($merge_contig->name, \@merge_contig_tags);
    $merge_contig->tags(\@merge_contig_tags	);

    return $merge_contig;
}




sub _get_ace_object
{
    Carp::confess("This method (_get_ace_object) was removed because it interacted with GSC classes. Contact $ENV{GENOME_EMAIL_PIPELINE} if you need this functionality.");
}

sub _get_phd_object
{
    Carp::confess("This method (_get_phd_object) was removed because it interacted with GSC classes. Contact $ENV{GENOME_EMAIL_PIPELINE} if you need this functionality.");
}

sub _calculate_split_region
{
    my ($left_reads, $right_reads, $contig_length) = @_;
    my $rRightMostLeftRead;
    my $rLeftMostLeftRead;
    my $nLeftMostLeftRead;
    my $nRightMostLeftRead;

    my $init = 0;
    foreach my $read (values %{$left_reads})
    {
        if ( $init == 0 )
        {
            $nRightMostLeftRead = $read->align_end;
            $rRightMostLeftRead = $read;
            $nLeftMostLeftRead  = $read->align_start;
            $rLeftMostLeftRead = $read;
            $init = 1;
        }
        else
        {
            if ( $nRightMostLeftRead < $read->align_end )
            {
                $nRightMostLeftRead = $read->align_end;
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

    foreach my $read (values %{$right_reads})
    {
        if ( $init == 0 )
        {
            $nLeftMostRightRead = $read->align_start;
            $rLeftMostRightRead  = $read;
            $nRightMostRightRead = $read->align_end;
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
            if ( $nRightMostRightRead < $read->align_end )
            {
                $nRightMostRightRead = $read->align_end;
                $rRightMostRightRead = $read;
            }
        }
    }

    # now find the overlap region
    #
    # the assertion is due to the fact that, since the assembly was
    # original contiguous, there is no way that it can be divided into
    # two groups of reads that don't overlap

    # the start_pos_, nRightEndOfTear_ region marks the torn region--the
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
    #			 start_pos_	   nRightEndOfTear_
    # It is NOT:
    #		^											 ^
    #		here			  to						 here
    # That is, a read that is at the cursor where the user is tearing,
    # may not be completely contained within the torn region.

    my ($start_pos, $end_pos);
    assert(
        calculate_intersection(
            $nLeftMostRightRead, $nRightMostRightRead,
            $nLeftMostLeftRead,  $nRightMostLeftRead,
            \$start_pos,    \$end_pos
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
    $start_pos = 1 if($start_pos < 1);
    $end_pos = $contig_length if($end_pos > $contig_length);
    return ($start_pos, $end_pos);
}

sub get_user_selected_reads
{
    my ($ref_reads) = @_;
    init_gtk() unless $initialized;
    Gtk->set_locale; 
    Gtk->init;

    my %reads = %{$ref_reads};

    my $d = Gtk::Dialog->new();
    $d->set_usize(300,300);
    $d->show;

    my $sw = Gtk::ScrolledWindow->new();
    $sw->set_policy("automatic", "automatic");
    #$sw->set_policy("always", "always");
    $sw->show;
    $d->vbox->pack_start_defaults($sw);

    my $vbox = Gtk::VBox->new();
    $vbox->show;
    $sw->add_with_viewport($vbox);

    foreach my $read_name (keys %reads)
    {
        my $bc = Gtk::ButtonCrate::Radio->new('h');
        $vbox->pack_start($bc->bbox, 0, 0, 0);

        my $label = Gtk::Label->new($read_name);
        $label->show;
        $bc->bbox->pack_start($label, 0, 0, 0);

        $bc->add_buttons(" <-", " ->");

        $bc->set_button_active_by_num(1) if($reads{$read_name} eq "right");	

        $bc->bbox->set_usize(100, 50);

        $reads{$read_name} = $bc;
    }

    my $b = Gtk::Button->new("Done");
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
            Gtk->main_quit;
        },
        \%reads,
        $ref_reads,
        $d
    );
    $b->show;
    $d->action_area->pack_start_defaults($b);

    $b = Gtk::Button->new("Cancel");
    $b->signal_connect("clicked", sub{ Gtk->exit(0) });
    $b->show;
    $d->action_area->pack_start_defaults($b);

    Gtk->main;

}
sub split
{
    my ($self, $contig, $phd_object, %params) = @_;		

    my $split_position = $params{split_position} if exists $params{split_position};

    my $in_fh;
    my $out_fh;

    my $bKeepCounting=1;
    my $no_gui=0;
    $no_gui = $params{'no_gui'} if exists $params{'no_gui'};		

    my %user_select_reads;
    my $padded_split_position = $contig->_transform->get_pad_position($split_position);
    my %reads = %{$contig->reads};
    my %left_reads;
    my %right_reads;
    foreach my $read (values %reads) 
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

        my $read_position = ($read->get_parent_position_from_child_position ($read->qual_clip_start, 1) + $read->get_parent_position_from_child_position($read->qual_clip_end, 1))/2.0;
        if($no_gui)
        {
            if ( $read_position < $padded_split_position )
            {
                $left_reads{$read->name} = $read;
            }
            else
            {
                $right_reads{$read->name} = $read;
            }
        }
        else
        {
            if ( $read->align_end < $padded_split_position)
            {
                $left_reads{$read->name} = $read;
            }
            elsif( $read->align_start > $padded_split_position )
            {
                $right_reads{$read->name} = $read;
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
        %left_reads = %{$params{left_reads}};
        %right_reads = %{$params{right_reads}};
    }

    if(keys %user_select_reads)
    {
        get_user_selected_reads(\%user_select_reads);
        foreach my $read_name (keys %user_select_reads)
        {
            if($user_select_reads{$read_name} eq "left")
            {
                $left_reads{$read_name} = $reads{$read_name};
            }
            else
            {
                $right_reads{$read_name} = $reads{$read_name};
            }	
        }
    }
    assert( ( scalar keys (%left_reads) + scalar keys (%right_reads) ) == scalar keys %reads );	
    die(   "You are putting all of the reads into the right contig so this is not a tear") if ( (scalar keys %left_reads) == 0 ) ;
    die(   "You are putting all of the reads into the left contig so this is not a tear") if ( (scalar keys %right_reads) == 0 );

    # at this point, all reads in the contig have been put into 2 lists:
    # a list of reads that the user wants to go into the 'left' contig,
    # and a list of reads that the user wants to go into the 'right contig.
    # (Note that these contigs themselves may fall apart into more contigs.)

    # Find the torn region--To the left and right of this region, any column
    # will be precisely the same as it was before the tear:  the same reads
    #  be in that column.  Within the tear, the reads will be different,
    # and thus we will have to recompute the consensus, Golden path, and
    # base segment.

    my ($start_pos, $end_pos) = _calculate_split_region(\%left_reads,\%right_reads,$contig->length);
    print "tearing contigs from ", $contig->_transform->get_unpad_position($start_pos), 
    " to ", $contig->_transform->get_unpad_position($end_pos), "\n";


    # now lets find if things fall into more than 2 contigs in the torn region
    # To do this, let's make two little lists:  The left reads in the
    # torn region, and the right reads in the torn region.
    my $left_overlap_reads = $contig->get_overlapping_reads($start_pos, $end_pos, \%left_reads);
    my $right_overlap_reads = $contig->get_overlapping_reads($start_pos, $end_pos, \%right_reads);

    foreach( values %{$left_overlap_reads},values %{$right_overlap_reads})
    {
        my $phd = $phd_object->get_phd($_->phd_file);
        my @quality = @{$phd->unpadded_base_quality};
        @quality = reverse @quality if ($_->complemented);
        $_->unpadded_base_quality( \@quality);				
    }
    # Now let's check (for the left contig and the right contig)
    # that there is a base for every position
    # (unaligned or aligned) for this entire region

    exit if ( !check_contiguous( [values %{$left_overlap_reads}], $start_pos, $end_pos ) || 
        !check_contiguous( [values %{$right_overlap_reads}], $start_pos, $end_pos ));	

    # point of no return
    # At this point, there is no longer any error checking.
    # Errors indicate a program bug.

    # Transfer reads to new contigs.
    # Fix alignment positions and tag positions of the right contig 326-592, 592-870

    my $left_contig = Genome::Model::Tools::Pcap::Contig->new;
    $left_contig->name( $contig->name . "a" );
    $left_contig->reads( \%left_reads);
    $left_contig->complemented( $contig->complemented);

    # in the torn region, figure out the consensus base and quality
    # The fast method of doing this is to go through each read once,
    # and mark two arrays--the highest quality at a cons pos and the
    # read that attained that highest quality.  I will ignore the
    # "aligned" clipping at this point, since a read could be high
    # quality unaligned precisely because there is a mis-join, which
    # is the reason the user is tearing the contig.  Thus the correct
    # consensus could be this unaligned high quality region.

    my $left_cont_start_pos = max( $start_pos,  1 );
    my $left_cont_end_pos   = min( $end_pos, $contig->base_count );
    assert( $left_cont_start_pos <= $left_cont_end_pos );
    $left_contig->padded_base_string(substr( $contig->padded_base_string, 0, $left_cont_end_pos ));
    my @temp_qual = @{$contig->padded_base_quality};
    $left_contig->padded_base_quality([ splice(@temp_qual, 0, $left_cont_end_pos)]);
    my ($best_quality_reads, $best_quality) = $left_contig->recalculate_consensus($start_pos,$end_pos,$left_overlap_reads);

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
    if ( 1 < $start_pos )
    {
        my $nAfterLastBaseSegmentToBeTransferred =
        $contig->get_bs_index_at_position( $start_pos );

        if ( $nAfterLastBaseSegmentToBeTransferred == -1 )
        {
            print "nAfterLastBaseSegmentToBeTransferred == -1, start_pos = ",
            $start_pos, "\n";
            assert(0);
        }

        my $nLastBaseSegmentToBeTransferred =
        $nAfterLastBaseSegmentToBeTransferred - 1;

        # note that if the base segment that covers start_pos_ is the
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
                start_pos => $rBaseSeg->{start_pos},
                end_pos   => $rBaseSeg->{end_pos},
                read_name => $rBaseSeg->{read_name},								
            };
            push( @left_bs, $rNewBaseSeg );
        }

        # now look at the base segment that perhaps must be cut short
        # We now guarantee that this base segment (or the one we just added),
        # will end precisely before the tear region

        my $rBoundaryBaseSeg = $bs[ $nLastBaseSegmentToBeTransferred + 1 ];

        if ( $rBoundaryBaseSeg->{start_pos} < $start_pos )
        {
            assert( $rBoundaryBaseSeg->{end_pos} >= $start_pos );
            my $rNewBoundaryBaseSeg = {
                type      => 'base_segment',
                start_pos => $rBoundaryBaseSeg->{start_pos},
                end_pos   => ( $start_pos - 1 ),
                read_name => $rBoundaryBaseSeg->{read_name},
            };

            push( @left_bs, $rNewBoundaryBaseSeg );
        }
    }

    # now add base segments for the tear region
    # We will add a new base segment for each base and not worry about
    # trying to consolidate adjacent base segments that share the same
    # read.

    # Fixed Bug:  DG, October 4, 2000.  It is possible that this torn region
    # can extend beyond the end of the consensus.  Clearly it doesn't
    # make sense to have base segments there, in either new contig.  So
    # only create base segments in the torn region that overlaps the consensus.

    my $nTornRegionOnConsensusLeft;
    my $nTornRegionOnConsensusRight;

    if (
        calculate_intersection(
            $start_pos,
            $end_pos,
            1,                              #start of consensus
            $left_contig->base_count,         #end of consensus
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
                start_pos => $pos,
                end_pos   => $pos,
                read_name => $best_quality_reads->[$pos]->name,
            };
            push( @left_bs, $rNewBaseSegment );
        }
    }
    $left_contig->base_segments( \@left_bs);

    # now recalculate the consensus quality values for the new left contig
    # (The reason for the intersect is that there may be reads that
    # stick out on the left to negative consensus positions and some
    # of these may be in the right contig and thus the start_pos_
    # may be to the left of where the new left contig starts.)

    my $nConsPosLeftToRecalculate;
    my $nConsPosRightToRecalculate;

    if (
        calculate_intersection(
            1,                      $left_contig->base_count,
            $start_pos,             $end_pos,
            \$nConsPosLeftToRecalculate, \$nConsPosRightToRecalculate
        )
    )
    {
        Genome::Model::Tools::Pcap::Utility::recalculate_consensus_qualities_and_change($left_contig, $nConsPosLeftToRecalculate, $nConsPosRightToRecalculate, [values %{$left_overlap_reads}]);
    }

    # That completes the base segments for the new left contig

    # Now work on the new right contig.

    # add reads to new contig

    my $offset = -$start_pos + 1;

    # handle case like this:
    #			cccccccccccccccccccccccccc
    #	---------------------------
    #		 +++++++++++++++++++++++++++++++++++

    # where the --- go into the new left contig and the
    # +++ go into the new right contig.
    # In this case, the beginning of the ccc... is 1 in
    # both the new left contig and the new right contig
    # Thus there is no translation.

    if ( $start_pos < 1 )
    {
        $offset = 0;
    }

    #now we recalculate the coordinates for the reads that are going into the right
    #contig
    foreach my $read (values %right_reads)
    {
        $read->position($read->position + $offset);
    }

    my $right_contig = Genome::Model::Tools::Pcap::Contig->new;
    $right_contig->name( $contig->name . "b" );
    $right_contig->reads( \%right_reads);
    $right_contig->complemented( $contig->complemented);

    # reset the consensus base and quality arrays to refigure them
    # based on the reads from the right contig.  Note that it is
    # possible (likely!) that the consensus bases and quality in the
    # torn region will be different for the two ends of the new contigs

    my $right_cont_start_pos = max( $start_pos,  1 );
    my $right_cont_end_pos = min( $end_pos, $contig->base_count );
    assert( $right_cont_start_pos <= $right_cont_end_pos );
    $right_contig->padded_base_string(substr( $contig->padded_base_string, $right_cont_start_pos-1, $contig->base_count-$right_cont_start_pos + 1 ));
    @temp_qual = @{$contig->padded_base_quality};
    $right_contig->padded_base_quality([ splice(@temp_qual, $right_cont_start_pos-1, $contig->base_count-$right_cont_start_pos + 1)]);

    $right_cont_start_pos += $offset;
    $right_cont_end_pos += $offset;
    ($best_quality_reads, $best_quality) = $right_contig->recalculate_consensus($right_cont_start_pos,$right_cont_end_pos,$right_overlap_reads);

    # now create base segments for the torn region

    # Fixed Bug:  DG, October 4, 2000.  It is possible that this torn region
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

    if ( $right_cont_start_pos < $right_cont_end_pos )
    {
        @right_bs = @{$right_contig->calculate_base_segments($right_cont_start_pos, $right_cont_end_pos,$best_quality_reads)};
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

    if ( $end_pos < $contig->base_count )
    {
        my $nBoundaryBaseSegment =
        $contig->get_bs_index_at_position($end_pos +1, \@bs,  );

        my $rBoundaryBaseSegment = $bs[$nBoundaryBaseSegment];

        # create a new base segment which starts precisely at
        # end_pos + 1.  This may be the same as this
        # base segment, or it may be to the right of it.
        my $rNewBoundaryBaseSegment = {
            type      => 'base_segment',
            start_pos => $end_pos + 1 + $offset,
            end_pos => ($rBoundaryBaseSegment->{end_pos} + $offset),
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
                start_pos => $rOldBaseSeg->{start_pos} + $offset,
                end_pos => ( $rOldBaseSeg->{end_pos} + $offset ),
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

    my $clipped_start_pos;
    my $clipped_end_pos;
    if (
        calculate_intersection(
            $right_cont_start_pos,
            $right_cont_end_pos,
            1,
            $right_contig->base_count,
            \$clipped_start_pos,
            \$clipped_end_pos
        )
    )
    {
        Genome::Model::Tools::Pcap::Utility::recalculate_consensus_qualities_and_change($right_contig,
            $clipped_start_pos,
            $clipped_end_pos,
            [values %right_reads]);
    }

    _transfer_consensus_tags($contig,$left_contig,$right_contig,$start_pos,$end_pos, $offset);

    printf("tear complete.\n");    	

    return ($left_contig,$right_contig);		
}

sub _transfer_consensus_tags
{
    my ($contig, $left_contig, $right_contig, $start_pos, $end_pos, $offset) = @_;
    my @left_tags,
    my @right_tags;
    foreach my $tag (@{$contig->tags})
    {

        if (    $tag->stop <= $start_pos
            && $tag->start < $end_pos )
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
        elsif (    $start_pos <= $tag->start
            && $end_pos < $tag->stop )
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

            my $right_tag = $tag->copy($tag);			

            my $left_tag = $tag;

            my $nOldConsPosEnd = $tag->stop;
            $left_tag->stop( min( $end_pos, $tag->stop ));
            $left_tag->parent($left_tag->parent .'a');

            $right_tag->start( max( $right_tag->start, $start_pos ));
            $right_tag->start($right_tag->start + $offset);

            $right_tag->stop( $nOldConsPosEnd + $offset);

            $right_tag->parent($right_tag->parent . 'b');

            printf("had to split tag of type %s into two:  one from %d to %d in contig %s and the other from %d to %d in contig %s",
                $tag->type,                 $left_tag->start,
                $left_tag->stop,    $left_contig->name,
                $right_tag->start, $right_tag->stop,
                $right_contig->name
            );
        }
    }
    $left_contig->tags(\@left_tags);
    $right_contig->tags(\@right_tags);
}


#This split functions handles all fileio, PSE creation, etc.
sub read_data_and_split_contig
{
    Carp::confess("This method (read_data_and_split_contig) was removed because it interacted with GSC classes. Contact $ENV{GENOME_EMAIL_PIPELINE} if you need this functionality.");
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
    $contig->unpadded_base_quality(reverse @{$contig->unpadded_base_quality});

    #next complement reads and read positions for each read in the read array
    my $reads = $contig->children;
    foreach my $read (values %$reads)
    {
        my $temp = reverse $read->padded_base_string;
        $read->padded_base_string ($temp =~ tr/actgACTG/tgacTGAC/);

        $read->position(  $nConsLength - ($read->position + $read->length)+2);
        $read->complemented(!$read->complemented);

        my @aReadTags = @{$read->tags};
        #complement read tags
        foreach my $tag (@aReadTags)
        {
            my $nReadLength = $read->length;
            #first change the offset, so that the start position and end
            #positions are now offsets from the end of the contig instead
            #of the beginning of the contig
            $tag->start_pos( $nReadLength - $tag->start_pos+1); 
            $tag->end_pos ($nReadLength - $tag->end_pos+1);
            #now swap start and end positions
            my $temp = $tag->end_pos;
            $tag->end_pos( $tag->start_pos);
            $tag->start_pos($temp); 	
        }
        $read->tags(\@aReadTags);
    }

    #now complement base segments
    my $base_segs = $contig->base_segments;
    foreach my $bs (@{$base_segs})
    {
        $bs->{start_pos} = $nConsLength - $bs->{start_pos}+1; 
        $bs->{end_pos} = $nConsLength - $bs->{end_pos}+1;
        my $temp = $bs->{end_pos};
        $bs->{end_pos} = $bs->{start_pos};
        $bs->{start_pos} = $temp;
    }

    @{$base_segs} = reverse @{$base_segs};
    $contig->base_segments($base_segs);



    #complement contig tags, if there are any for this particular contig
    my @aContigTags = @{$contig->tags};
    foreach my $tag (@aContigTags)
    {
        $tag->start_pos( $nConsLength - $tag->start_pos+1);
        $tag->end_pos ( $nConsLength - $tag->end_pos+1);

        #now swap start and end positions
        my $temp = $tag->end_pos;
        $tag->end_pos ( $tag->start_pos );
        $tag->start_pos ( $temp);	
    }
    $contig->tags(\@aContigTags);		
}

sub read_data_and_complement_contig
{
    Carp::confess("This method (read_data_and_complement_contig) was removed because it interacted with GSC classes. Contact $ENV{GENOME_EMAIL_PIPELINE} if you need this functionality.");
}

sub read_data_and_merge_contigs
{
    Carp::confess("This method (read_data_and_merge_contigs) was removed because it interacted with GSC classes. Contact $ENV{GENOME_EMAIL_PIPELINE} if you need this functionality.");
}

sub remove_reads
{
    my ($self, $contig, $phd_object, %params) = @_;		

    my %reads = %{$contig->reads};
    my $remove_reads = $params{remove_reads};
    foreach my $read_name (@{$remove_reads}) 
    { 
        if(delete $reads{$read_name})
        {
            print "Successfully removed $read_name from ".$contig->name."\n";
        }
        else
        {
            warn "Could not find $read_name in ".$contig->name."\n";
        }
    }
    if((scalar keys %reads) == 0)
    {
        print $contig->name." has no reads left.\n";
        return undef;
    }	
    my $region = $contig->get_contiguous_reads(\%reads,$contig->tags);
    if(!defined $region) {print $contig->name."\n"; return undef;}
    $contig->reads($region->{reads});%reads = %{$region->{reads}};
    $contig->tags($region->{tags});
    # at this point, all reads in the contig have been put into 2 lists:
    # a list of reads that the user wants to go into the 'left' contig,
    # and a list of reads that the user wants to go into the 'right contig.
    # (Note that these contigs themselves may fall apart into more contigs.)

    # Find the torn region--To the left and right of this region, any column
    # will be precisely the same as it was before the tear:  the same reads
    #  be in that column.  Within the tear, the reads will be different,
    # and thus we will have to recompute the consensus, Golden path, and
    # base segment.

    my $start_pos = 1;
    my $end_pos = $region->{length};
    foreach( values %reads)
    {
    	my $phd = $phd_object->get_phd($_->phd_file);
    	my @quality = @{$phd->unpadded_base_quality};
        @quality = reverse @quality if ($_->complemented);
    	$_->unpadded_base_quality( \@quality);
    }

    print "Done loading read quality for ".$contig->name."\n";
    exit if ( !check_contiguous( [values %reads], $start_pos, $end_pos ) );	

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


    assert( $start_pos <= $end_pos );
    my ($best_quality_reads, $best_quality) = $contig->calculate_consensus($start_pos,$end_pos,\%reads);

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

    $contig->base_segments($contig->calculate_base_segments($start_pos, $end_pos,$best_quality_reads));

    # now recalculate the consensus quality values for the new left contig
    # (The reason for the intersect is that there may be reads that
    # stick out on the left to negative consensus positions and some
    # of these may be in the right contig and thus the start_pos_
    # may be to the left of where the new left contig starts.)

    Genome::Model::Tools::Pcap::Utility::recalculate_consensus_qualities_and_change($contig, $start_pos, $end_pos, [values %reads]);

    # That completes the base segments for the new left contig
    # do full checking once
    assert( $contig->is_bs_array_structure_ok( 1 ));

    printf("read removal complete.\n");    	

    return $contig;		
}

sub remove_read
{
    my ($self, $contig, $phd_object, %params) = @_;		

    my %user_select_reads;
    my %reads = %{$contig->reads};
    my $read = delete $reads{$params{remove_read}};
    $contig->reads(\%reads);
    # at this point, all reads in the contig have been put into 2 lists:
    # a list of reads that the user wants to go into the 'left' contig,
    # and a list of reads that the user wants to go into the 'right contig.
    # (Note that these contigs themselves may fall apart into more contigs.)

    # Find the torn region--To the left and right of this region, any column
    # will be precisely the same as it was before the tear:  the same reads
    #  be in that column.  Within the tear, the reads will be different,
    # and thus we will have to recompute the consensus, Golden path, and
    # base segment.

    my $start_pos = $read->align_start;
    my $end_pos = $read->align_end;	
    # now lets find if things fall into more than 2 contigs in the torn region
    # To do this, let's make two little lists:  The left reads in the
    # torn region, and the right reads in the torn region.
    my $overlap_reads = $contig->get_overlapping_reads($start_pos, $end_pos, \%reads);

    foreach( values %{$overlap_reads})
    {
        my $phd = $phd_object->get_phd($_->phd_file);
        my @quality = @{$phd->unpadded_base_quality};
        @quality = reverse @quality if ($_->complemented);
        $_->unpadded_base_quality( \@quality);				
    }
    # Now let's check (for the left contig and the right contig)
    # that there is a base for every position
    # (unaligned or aligned) for this entire region
    my $cont_start_pos = max( $start_pos,  1 );
    my $cont_end_pos   = min( $end_pos, $contig->base_count );
    exit if ( !check_contiguous( [values %{$overlap_reads}], $cont_start_pos, $cont_end_pos ) );	

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


    assert( $cont_start_pos <= $cont_end_pos );
    my ($best_quality_reads, $best_quality) = $contig->recalculate_consensus($cont_start_pos,$cont_end_pos,$overlap_reads);

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

    my @bs = @{$contig->base_segments};
    my @bs_new = @{$contig->calculate_base_segments($cont_start_pos, $cont_end_pos,$best_quality_reads)};
    my $bs_start = $contig->get_bs_index_at_position($cont_start_pos);
    my $bs_end = $contig->get_bs_index_at_position($cont_end_pos);
    $bs[$bs_start-1]->{end_pos} = $cont_start_pos-1 if($cont_start_pos>1);
    $bs[$bs_end+1]->{start_pos} = $cont_end_pos+1 if($cont_end_pos<$contig->base_count);
    splice(@bs, $bs_start,$bs_end-$bs_start+1,@bs_new);
    $contig->base_segments( \@bs);

    # now recalculate the consensus quality values for the new left contig
    # (The reason for the intersect is that there may be reads that
    # stick out on the left to negative consensus positions and some
    # of these may be in the right contig and thus the start_pos_
    # may be to the left of where the new left contig starts.)

    my $nConsPosLeftToRecalculate;
    my $nConsPosRightToRecalculate;

    if (
        calculate_intersection(
            1,                      $contig->base_count,
            $start_pos,             $end_pos,
            \$nConsPosLeftToRecalculate, \$nConsPosRightToRecalculate
        )
    )
    {
        Genome::Model::Tools::Pcap::Utility::recalculate_consensus_qualities_and_change($contig, $nConsPosLeftToRecalculate, $nConsPosRightToRecalculate, [values %{$overlap_reads}]);
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

    print "No x/ns found in contig: " . $contig->name . "\n" and return undef
    unless scalar @{$xn_positions} > 0;

    #GET THE ACTUAL BASE AT THAT POSITION FROM EACH READ
    #ALIGNED UNDER THE N/X POSITION

    my $reads = $contig->reads;
    my $base_counts = {};
    my $replace_xns = {};

    foreach my $pos (@$xn_positions)
    {
        foreach my $read (values %{$reads})
        {
            my @read_positions = ($read->start_position .. $read->end_position);
            next unless grep (/^$pos$/, @read_positions);
            my $read_position = $pos - $read->start_position;
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

Genome::Model::Tools::Pcap::ContigTools performs operation on objects of type contig.  So far the functions supported are merging and splitting.

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
