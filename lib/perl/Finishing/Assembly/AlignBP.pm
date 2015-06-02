package Finishing::Assembly::AlignBP;
our $VERSION = 0.01;

use Finishing::Assembly::Transform;
use Carp::Assert;
use Bio::Seq;
use Bio::Tools::dpAlign;
use strict;
use warnings;

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
	my ($left_string, $right_string, %params) = @_;
	
	my $gap_open_penalty = $params{gap_open_penalty};
	my $gap_ext_penalty = $params{gap_ext_penalty};
	my $max_length = $params{max_length};
	my $left_align = { original_string => $left_string };
	my $right_align = { original_string => $right_string };
	$left_string =~ tr/*/r/;
	$right_string =~ tr/*/r/;
	if(length($left_string) >= $max_length)
	{
		$left_string = substr($left_string, length($left_string)-$max_length, $max_length);	
	}

	if(length($right_string) >= $max_length)
	{
		$right_string = substr($right_string, 0, $max_length);	
	}

	my $left_seq = Bio::Seq->new( -display_id => "left_string",
		-seq => $left_string

	);

	my $right_seq = Bio::Seq->new( -display_id => "right_string",
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
	_convert_bioperl_to_gsc_transform($out,$left_align,$right_align,%params);
	return ($left_align,$right_align);

}
#align object
#	start - start of alignment in orignal string coordinates
#   end - end of alignment in original string coordinates
#	local_align_seq - locally aligned sequence
#	align_seq - entire sequence containing newly aligned region
#	merge_seq - sequence containing part of string1 to left of align region, and part of string2 to the right of align region
#	transform


sub _convert_bioperl_to_gsc_transform
{
	my($out,$left_align,$right_align,%params) = @_;
	my $max_length = $params{max_length};
	#we'll go ahead and create a hash named 'align' which will contain all of the
	#data that we are going to use when performing our alignment
	#think of it as a convenient way of passing around the alignment data for each
	#contig
	$left_align->{start}= $out->{_start_end_lists}{"left_string"}[0]{start};
	$left_align->{end} = $out->{_start_end_lists}{"left_string"}[0]{end};

	#appears in the left contig to the one produced by bioperl
	#we use it later when adding these pads to the reads	
	$left_align->{local_align_seq} = $left_align->{align_seq} = $out->{_start_end_lists}{"left_string"}[0]{seq};
	$right_align->{local_align_seq} = $right_align->{align_seq} = $out->{_start_end_lists}{"right_string"}[0]{seq};
	$left_align->{align_seq} =~ tr/rR/**/;
	$right_align->{align_seq} =~ tr/rR/**/;
	if( length ($left_align->{original_string}) > $max_length)
	{
		$left_align->{start} = length ($left_align->{original_string})-$max_length+$left_align->{start};
		$left_align->{end} = length ($left_align->{original_string})-$max_length+$left_align->{end};
	}
	$right_align->{start}= $out->{_start_end_lists}{"right_string"}[0]{start};
	$right_align->{end} = $out->{_start_end_lists}{"right_string"}[0]{end};
	my $left_pre = substr( $left_align->{original_string}, 0, $left_align->{start}-1 );
	my $left_post = substr( $left_align->{original_string}, $left_align->{end}, length ($left_align->{original_string})-$left_align->{end});
	$left_align->{merge_seq} = $left_pre.$left_align->{align_seq};
	$left_align->{align_seq} = $left_pre.$left_align->{align_seq}.$left_post;	
	my $right_pre = substr( $right_align->{original_string}, 0, $right_align->{start}-1 );
	my $right_post = substr( $right_align->{original_string}, $right_align->{end}, length ($right_align->{original_string})-$right_align->{end});
	$right_align->{align_seq} = $right_pre.$right_align->{align_seq}.$right_post;
	$left_align->{merge_seq} .= $right_post;
	$left_align->{merge_seq} =~ tr/-/*/;


	my $left_transform = Finishing::Assembly::Transform->new($left_align->{align_seq}, '-');#this creates a transform for the merge region as it originally 
	my $right_transform = Finishing::Assembly::Transform->new($right_align->{align_seq}, '-');#this creates a transform for the merge region as it originally 
	$right_transform->_offset( $left_align->{start}-$right_align->{start});
	#may want to check if offset is valid
	$left_align->{transform} = $left_transform;
	$right_align->{transform} = $right_transform;
}

1;
