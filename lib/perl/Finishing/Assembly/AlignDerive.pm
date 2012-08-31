package Finishing::Assembly::AlignDerive;
our $VERSION = 0.01;

use Finishing::Assembly::Transform;
use Carp::Assert;
use strict;
use warnings;
use List::Util qw(min max);

#sub new 
#{
#    croak("__PACKAGE__:new:no class given, quitting") if @_ < 1;
#	my ($caller, %params) = @_; 
#    my $caller_is_obj = ref($caller);
#    my $class = $caller_is_obj || $caller;
#    my $self = {};
#    bless ($self, $class);    
#    return $self;
#}

sub align
{
	my ($left_contig, $right_contig, %params) = @_;
	my $read_name = $params{read_name};
	my $left_reads = $left_contig->reads;
	my $right_reads = $right_contig->reads;
	my $left_read = $left_reads->{$read_name};
	my $right_read = $right_reads->{$read_name};
	my $left_transform = $left_read->sequence->get_transform;
	my $right_transform = $right_read->sequence->get_transform;
	my $new_left_transform = Finishing::Assembly::Transform->new;
	$new_left_transform->derive_from_transforms($left_transform,$right_transform);
	my $new_right_transform = Finishing::Assembly::Transform->new;
	$new_right_transform->derive_from_transforms($right_transform,$left_transform);
	$params{left_read} = $left_read;
	$params{right_read} = $right_read;
	$params{new_left_transform} = $new_left_transform;
	$params{new_right_transform} = $new_right_transform;
	$params{left_contig} = $left_contig;
	$params{right_contig} = $right_contig;
	my $left_align = { original_string => $left_contig->sequence->padded_base_string };
	my $right_align = { original_string => $right_contig->sequence->padded_base_string };
	_create_align_objects($left_align, $right_align, \%params);
	return ($left_align, $right_align);
}

sub _create_align_objects
{
	my ($left_align,$right_align,$params) = @_;
	my $left_read = $params->{left_read};
	my $right_read = $params->{right_read};
	my $left_read_length = $left_read->length;
	my $right_read_length = $right_read->length;
	my $new_left_transform = $params->{new_left_transform};
	my $new_right_transform = $params->{new_right_transform};
	my $new_left_read_length = $new_left_transform->pad_length;
	my $new_right_read_length = $new_right_transform->pad_length;
	
	#we'll go ahead and create a hash named 'align' which will contain all of the
	#data that we are going to use when performing our alignment
	#think of it as a convenient way of passing around the alignment data for each
	#contig
	$left_align->{start}= $left_read->position;
	$left_align->{end} = $left_read->end_position;

	$right_align->{start}= $right_read->position;
	$right_align->{end} = $right_read->end_position;
	
	($left_align->{align_seq},$right_align->{align_seq}) = get_align_seq2($left_align, $right_align, $params);
	#clip_align_seq($left_align, $right_align, $params);	
	$left_align->{local_align_seq} = $left_align->{align_seq};
	$right_align->{local_align_seq} = $right_align->{align_seq};
	#$left_align->{local_align_seq} = $left_align->{align_seq} = $new_left_transform->pad_string($left_read->sequence->padded_base_string);
	#$right_align->{local_align_seq} = $right_align->{align_seq} = $new_right_transform->pad_string($right_read->sequence->padded_base_string);
	
	my $left_pre = substr( $left_align->{original_string}, 0, $left_align->{start}-1 );
	my $left_post = substr( $left_align->{original_string}, $left_align->{end}, length ($left_align->{original_string})-$left_align->{end});
	$left_align->{merge_seq} = $left_pre.$left_align->{align_seq};
	$left_align->{align_seq} = $left_pre.$left_align->{align_seq}.$left_post;	
	my $right_pre = substr( $right_align->{original_string}, 0, $right_align->{start}-1 );
	my $right_post = substr( $right_align->{original_string}, $right_align->{end}, length ($right_align->{original_string})-$right_align->{end});
	$right_align->{align_seq} = $right_pre.$right_align->{align_seq}.$right_post;
	$left_align->{merge_seq} .= $right_post;
	$left_align->{merge_seq} =~ tr/-/*/;
	#print "Original left consensus: $left_align->{original_string}\n";
	#print "Original right consensus: $right_align->{original_string}\n";
	#print "Left align seq: $left_align->{align_seq}\n";
	#print "Right align seq: $right_align->{align_seq}\n";
	#rint "Merge seq: $left_align->{merge_seq}\n";

	my $left_transform = Finishing::Assembly::Transform->new($left_align->{align_seq}, '-');#this creates a transform for the merge region as it originally 
	my $right_transform = Finishing::Assembly::Transform->new($right_align->{align_seq}, '-');#this creates a transform for the merge region as it originally 
	$right_transform->_offset( $left_align->{start}-$right_align->{start});
	#$left_align->{end} = $params->{new_left_transform}->get_pad_position($left_align->{end}-$left_align->{start})+$left_align->{start};
	#$right_align->{end} = $params->{new_right_transform}->get_pad_position($right_align->{end}-$right_align->{start})+$right_align->{start};
	#may want to check if offset is valid
	$left_align->{transform} = $left_transform;
	$right_align->{transform} = $right_transform;
}
#need function that will deal with the fact that some of the read may be clipped off
#it will also need to get the consensus sequence that corresponds to the read sequence

sub get_align_seq
{
	my ($left_align, $right_align, $params) = @_;
	my $left_read = $params->{left_read};
	my $right_read = $params->{right_read};
	my $left_contig = $params->{left_contig};
	my $right_contig = $params->{right_contig};
	my $new_left_transform = $params->{new_left_transform};
	my $new_right_transform = $params->{new_right_transform};
	my $right_added_pads = max(0, 1 - $right_read->start_position);
	my $left_added_pads = max($left_read->end_position-$left_contig->length,0);
	$params->{right_added_pads} = $right_added_pads;
	$params->{left_added_pads} = $left_added_pads;
	my $left_consensus = $left_contig->sequence->padded_base_string;
	my $right_consensus =  $right_contig->sequence->padded_base_string;
	my $pad_total = $left_added_pads+$right_added_pads;
	$left_consensus = substr($left_consensus,$right_added_pads+$left_read->start_position-1,length($left_consensus)-$right_added_pads-$left_read->start_position-$left_added_pads+1);
	$right_consensus = substr($right_consensus,0,$right_read->end_position-$left_added_pads);	

	my $rsub_transform = $new_right_transform->get_sub_transform(unpad_start => $right_added_pads,unpad_end => $new_right_transform->unpad_length-$left_added_pads-1); 
	my $lsub_transform = $new_left_transform->get_sub_transform(unpad_start => $right_added_pads, unpad_end => $new_left_transform->unpad_length-$left_added_pads-1);

	my $left_align_seq = $lsub_transform->pad_string($left_consensus);
	my $right_align_seq = $rsub_transform->pad_string($right_consensus);
	
	return ($left_align_seq, $right_align_seq);
}

sub clip_align_seq
{
	my ($left_align, $right_align, $params) = @_;
	my $new_left_transform = $params->{new_left_transform};
	my $new_right_transform = $params->{new_right_transform};
	my $right_added_pads = $params->{right_added_pads};
	my $left_added_pads = $params->{left_added_pads};
	my $start_clip = max($new_left_transform->get_pad_position($params->{left_read}->align_clip_start),
					 $new_right_transform->get_pad_position($params->{right_read}->align_clip_start),
					 $right_added_pads+1)-1;
	my $end_clip = min($new_left_transform->pad_length -
					$new_left_transform->get_pad_position($params->{left_read}->align_clip_end),
					$new_right_transform->pad_length -	
					$new_right_transform->get_pad_position($params->{right_read}->align_clip_end),
					$left_added_pads+1);
	$start_clip -= $right_added_pads;
	$end_clip -= $left_added_pads;				
	$left_align->{start} += $start_clip; 
	$left_align->{end} -= $end_clip;
	$right_align->{start} += $start_clip;
	$right_align->{end} -= $end_clip;
	my $temp_left_end = $new_left_transform->pad_length - $end_clip -1;
	my $temp_right_end = $new_right_transform->pad_length - $end_clip -1;
	$params->{new_left_transform} = $new_left_transform->get_sub_transform(pad_start => $start_clip, pad_end=> $temp_left_end);
	$params->{new_right_transform} = $new_right_transform->get_sub_transform(pad_start => $start_clip, pad_end => $temp_right_end);
	$left_align->{align_seq} = substr($left_align->{align_seq},$start_clip,$temp_left_end-$start_clip+1);
	$right_align->{align_seq} = substr($right_align->{align_seq},$start_clip,$temp_right_end-$start_clip+1);
}

sub get_align_seq2
{
	my ($left_align, $right_align, $params) = @_;
	my $left_read = $params->{left_read};
	my $right_read = $params->{right_read};
	my $left_contig = $params->{left_contig};
	my $right_contig = $params->{right_contig};
	my $new_left_transform = $params->{new_left_transform};
	my $new_right_transform = $params->{new_right_transform};
	my $start_clip = max($new_right_transform->get_pad_position(max(1 - $right_read->start_position,0)),
						$new_left_transform->get_pad_position($params->{left_read}->align_clip_start)-1,
						$new_right_transform->get_pad_position($params->{right_read}->align_clip_start)-1);
	my $end_clip = max($new_left_transform->get_pad_position(max($left_read->end_position-$left_contig->length),0),
					$new_left_transform->pad_length -
					$new_left_transform->get_pad_position($params->{left_read}->align_clip_end),
					$new_right_transform->pad_length -	
					$new_right_transform->get_pad_position($params->{right_read}->align_clip_end));
	#unpadded start and end clip, since they're unpadded, they are specific to each read
	my $left_start_clip = $new_left_transform->get_unpad_position($start_clip);
	my $left_end_clip = $new_left_transform->unpad_length-$new_left_transform->get_unpad_position($new_left_transform->pad_length-$end_clip);
	my $right_start_clip = $new_right_transform->get_unpad_position($start_clip);
	my $right_end_clip = $new_right_transform->unpad_length-$new_right_transform->get_unpad_position($new_right_transform->pad_length-$end_clip);
	#					
	$left_align->{start} += $left_start_clip; 
	$left_align->{end} -= $left_end_clip;
	$right_align->{start} += $right_start_clip;
	$right_align->{end} -= $right_end_clip;
	my $left_consensus = $left_contig->sequence->padded_base_string;
	my $right_consensus =  $right_contig->sequence->padded_base_string;
	my $pad_total = $start_clip+$end_clip;
	my $clipped_length = $new_left_transform->pad_length - $pad_total;
	$left_consensus = substr($left_consensus,$left_align->{start}-1,$left_align->{end}-$left_align->{start}+1);
	$right_consensus = substr($right_consensus,$right_align->{start}-1,$right_align->{end}-$right_align->{start}+1);
	my $rsub_transform = $new_right_transform->get_sub_transform(unpad_start => $right_start_clip,unpad_end => $right_read->length-$right_end_clip-1); 
	my $lsub_transform = $new_left_transform->get_sub_transform(unpad_start => $left_start_clip, unpad_end => $left_read->length-$left_end_clip-1);
	#my $rsub_transform = $new_right_transform->get_sub_transform(pad_start => $start_clip,pad_end => $new_right_transform->pad_length-$end_clip); 
	#my $lsub_transform = $new_left_transform->get_sub_transform(pad_start => $start_clip, pad_end => $new_left_transform->pad_length-$end_clip);
	my $left_align_seq = $lsub_transform->pad_string($left_consensus,'-');
	my $right_align_seq = $rsub_transform->pad_string($right_consensus,'-');
	$params->{new_left_transform} = $lsub_transform;
	$params->{new_right_transform} = $rsub_transform;
	return ($left_align_seq, $right_align_seq);
}



1;
