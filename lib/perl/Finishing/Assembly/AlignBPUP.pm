package Finishing::Assembly::AlignBPUP;
our $VERSION = 0.01;

use Finishing::Assembly::Transform;
use Carp::Assert;
use Bio::Seq;
use Bio::Tools::dpAlign;
use List::Util qw(min max);
use strict;
use warnings;

sub align
{
	my ($left_string, $right_string, %params) = @_;
	
	my $gap_open_penalty = $params{gap_open_penalty};
	my $gap_ext_penalty = $params{gap_ext_penalty};
	my $max_length = $params{max_length};
	my $left_align = { original_string => $left_string };
	my $right_align = { original_string => $right_string };

	if(length($left_string) >= $max_length)
	{
		$left_string = substr($left_string, length($left_string)-$max_length, $max_length);	
	}

	if(length($right_string) >= $max_length)
	{
		$right_string = substr($right_string, 0, $max_length);	
	}
	
	my $left_transform = Finishing::Assembly::Transform->new($left_string);
	my $right_transform = Finishing::Assembly::Transform->new($right_string);
	my $left_string_padded = $left_string;
	my $right_string_padded = $right_string;
	$left_string =~ tr/*//d;
	$right_string =~ tr/*//d;
	my $left_seq = Bio::Seq->new( -display_id => "left_string",
		-seq => $left_string

	);

	my $right_seq = Bio::Seq->new( -display_id => "right_string",
		-seq => $right_string

	);

	$right_seq->alphabet('dna');
	$left_seq->alphabet('dna');
	my $factory = Bio::Tools::dpAlign->new ( -match => 1,
		-mismatch => -1,
		-gap => $gap_open_penalty,
		-ext => $gap_ext_penalty,
		-alg => Bio::Tools::dpAlign::DPALIGN_LOCAL_MILLER_MYERS);		
	
	my $out = $factory->pairwise_alignment($left_seq, $right_seq);
	my $left_seq_string = $out->{_start_end_lists}{"left_string"}[0]{seq};
	my $right_seq_string = $out->{_start_end_lists}{"right_string"}[0]{seq};
	my $bp_left_transform = Finishing::Assembly::Transform->new($left_seq_string,'-');
	my $bp_right_transform = Finishing::Assembly::Transform->new($right_seq_string,'-');
	my $clipped_left_transform = $left_transform->get_sub_transform(unpad_start => $out->{_start_end_lists}{"left_string"}[0]{start}-1,
																	unpad_end => $out->{_start_end_lists}{"left_string"}[0]{end}-1);
	my $clipped_right_transform = $right_transform->get_sub_transform(unpad_start => $out->{_start_end_lists}{"right_string"}[0]{start}-1,
																	  unpad_end => $out->{_start_end_lists}{"right_string"}[0]{end}-1);
	my $t1 = $clipped_left_transform;$t1->set_pad_char('a');
	my $t2 = $bp_left_transform;$t2->set_pad_char('b');
	my $t3 = $bp_right_transform;$t3->set_pad_char('c');
	my $t4 = $clipped_right_transform;$t4->set_pad_char('d');
	my $t5 = get_combined_transform($t1,$t2,$t3,$t4);
	my $t6 = get_combined_transform($t4,$t3,$t2,$t1);
	my $t7 = SIR(SIL($t5,$t2),$t1);
	my $t8 = SIR(SIL($t6,$t3),$t4);
	my $t9 = compact_pads($t7);
	my $t10	= compact_pads($t8);
	
	#need to convert output hash
	$out->{_start_end_lists}{"left_string"}[0]{start} = 
		$left_transform->get_pad_position($out->{_start_end_lists}{"left_string"}[0]{start}-1)+1;
	$out->{_start_end_lists}{"left_string"}[0]{end} =
		$left_transform->get_pad_position($out->{_start_end_lists}{"left_string"}[0]{end}-1)+1;
	$out->{_start_end_lists}{"right_string"}[0]{start} =
		$right_transform->get_pad_position($out->{_start_end_lists}{"right_string"}[0]{start}-1)+1;
	$out->{_start_end_lists}{"right_string"}[0]{end} =
		$right_transform->get_pad_position($out->{_start_end_lists}{"right_string"}[0]{end}-1)+1;
	
	$left_string_padded = substr
	(
		$left_string_padded,
		$out->{_start_end_lists}{"left_string"}[0]{start}-1,
		$out->{_start_end_lists}{"left_string"}[0]{end}-$out->{_start_end_lists}{"left_string"}[0]{start}+1
	);
	$right_string_padded = substr
	(
		$right_string_padded,
		$out->{_start_end_lists}{"right_string"}[0]{start}-1,
		$out->{_start_end_lists}{"right_string"}[0]{end}-$out->{_start_end_lists}{"right_string"}[0]{start}+1
	);
	$left_string_padded =~ tr/*/R/;
	$right_string_padded =~ tr/*/R/;
	my $new_left_seq = $t9->pad_string($left_string_padded,'-');
	my $new_right_seq = $t10->pad_string($right_string_padded,'-');
	
	$out->{_start_end_lists}{"left_string"}[0]{seq} = $new_left_seq;
	$out->{_start_end_lists}{"right_string"}[0]{seq} = $new_right_seq;
#	$out->{_start_end_lists}{"left_string"}[0]{seq} =~ tr/*/R/;
#	$out->{_start_end_lists}{"right_string"}[0]{seq} =~ tr/*/R/;
		
	_convert_bioperl_to_gsc_transform($out,$left_align,$right_align,%params);
	
	return ($left_align,$right_align);

}

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

sub get_combined_transform
{
	my ($t1, $t2, $t3, $t4) = @_;

	$t1 = $t1->copy($t1);
	$t2 = $t2->copy($t2);
	$t3 = $t3->copy($t3);
	$t4 = $t4->copy($t4);
	my $temp1 = combine_pads($t1,$t2);
	$temp1 = SIR($temp1,$t2);
	my $temp2 = combine_pads($t3,$t4);
	$temp2 = SIR($temp2,$t3);
	my $temp3 = combine_pads($temp1,$temp2);
	return $temp3;
}

#A1: [-----]
#A2: [---*--]
#B1: [----]
#B2: [--+--]
#C1: [----]
#C2: [--*--]
sub SL
{
	my ($t1, $t2) = @_;
	my $transform = Finishing::Assembly::Transform->new;
	my @t1ptu = @{$t1->{padded_to_unpadded}};
	my @t1utp = @{$t1->{unpadded_to_padded}};
	my @t2ptu = @{$t2->{padded_to_unpadded}};
	my @t2utp = @{$t2->{unpadded_to_padded}};
	
	#the easiest way to do this is to build a list of deletions
	my @del_list;
	my $pad_char = $t2->pad_char||'*';
	for(my $i = 0;$i<@t2ptu;$i++)
	{
		if($t2ptu[$i] eq $pad_char)
		{
			push @del_list, $i;
		}
	}
	@del_list = reverse sort @del_list;
	my $t1pad_char = $t1->pad_char;
	foreach(@del_list)
	{
		splice(@t1utp, $t1ptu[$_],1);
		for(my $i=$t1ptu[$_];$i<@t1utp;$i++){ $t1utp[$i]--; }
		splice(@t1ptu, $_,1);
		for(my $i=$_;$i<@t1ptu;$i++){ $t1ptu[$i]-- if($t1ptu[$i] =~/\d+/); }
	}
	
	$transform->{unpadded_to_padded} = \@t1utp;
	$transform->{padded_to_unpadded} = \@t1ptu;
	return $transform;

}

sub cmp_ins
{
	my ($a, $b) = @_;
	return 1 if $a->[0] > $b->[0];
	return -1 if $b->[0] > $a->[0];
	return 0;


}

sub SR
{
	my ($t1, $t2) = @_;
	my $transform = Finishing::Assembly::Transform->new;
	my @t1ptu = @{$t1->{padded_to_unpadded}};
	my @t1utp = @{$t1->{unpadded_to_padded}};
	my @t2ptu = @{$t2->{padded_to_unpadded}};
	my @t2utp = @{$t2->{unpadded_to_padded}};
	
	my $pad_char = $t2->pad_char;
	
	#the easiest way to do this is to build a list of insertions
	my @ins_list;
	for(my $i = 1;$i<@t2utp;$i++)
	{
		if(($t2utp[$i]-$t2utp[$i-1])>1)
		{
			push @ins_list, [$i-1, $t2utp[$i]-$t2utp[$i-1]-1];
		}
	}
	@ins_list = reverse sort { cmp_ins ($a, $b)} @ins_list;
	
	foreach(@ins_list)
	{
		#splice(@t1utp, $t1ptu[$_->[0]]+1,0,($_->[0]..($_->[0]+$_->[1]-1)) );
		#for(my $i=$t1ptu[$_->[0]];$i<@t1utp;$i++)
		#{ 
		#	$t1utp[$i]+=$_->[1] if($t1utp[$i] =~ /\d+/); 
		#}
		@t1utp = @t2utp;
		my $t1_position = $t1ptu[$_->[0]+$_->[1]];
		splice(@t1ptu, $_->[0]+1, 0, '+'x$_->[1] );#$t1utp[($_->[0]..($_->[0]+$_->[1]-1))]
		for(my $i=$_->[0];$i<@t1ptu;$i++)
		{ 
			if($t1ptu[$i] =~ /\d+/)
			{
				$t1ptu[$i]=$t1_position;
				$t1_position++;
			} 
		}
	}
	
	$transform->{unpadded_to_padded} = \@t1utp;
	$transform->{padded_to_unpadded} = \@t1ptu;
	return $transform; 
}

# a better description of what is happening is that the input side is having it's pads removed,
# and is then repadded with a bioperl string.
# These two operations only make sense when converting one kind of padded alignment to another.
#A1 : [acctg] [--*--]
#A2 : [acct+g] [--*-+-]
#B1 : [actg] [----]
#B2 : [ac*tg] [--*--]
#C1 : [ac*tg]  [----]
#C2 : [ac*t+g] [--*-+-]
#       -or-
#A1 : [acctg] [-----]
#A2 : [acct+g] [----+-]
#B1 : [actg] [----]
#B2 : [ac*tg] [--*--]
#C1 : [ac*tg]  [----]
#C2 : [ac*t+g] [--*-+-]

sub SIL
{
	my ($t1, $t2) = @_;
	my $transform = Finishing::Assembly::Transform->new;
	my @t1ptu = @{$t1->{padded_to_unpadded}};
	my @t1utp = @{$t1->{unpadded_to_padded}};
	my @t2ptu = @{$t2->{padded_to_unpadded}};
	my @t2utp = @{$t2->{unpadded_to_padded}};
	
	#the easiest way to do this is to build a list of deletions
	my @del_list;
	my $pad_char = $t2->pad_char;
	for(my $i = 0;$i<@t2ptu;$i++)
	{
		unless($t2ptu[$i] =~ /\d+/)
		{
			push @del_list, $i;
		}
	}
	@del_list = reverse sort @del_list;
	
	foreach(@del_list)
	{
		$t1ptu[$t1utp[$_]]=$pad_char;
		for(my $i=$t1utp[$_]+1;$i<@t1ptu;$i++){ $t1ptu[$i]-- if($t1ptu[$i] =~ /\d+/); }
		splice(@t1utp, $_,1);
		#for(my $i=$t1ptu[$_];$i<@t1utp;$i++){ $t1utp[$i]--; }
		
	}
	
	$transform->{unpadded_to_padded} = \@t1utp;
	$transform->{padded_to_unpadded} = \@t1ptu;
	return $transform;

}



#A1 : [acctg] [----]
#A2 : [acct+g] [--*-+-]
#B1 : [actg] [----]
#B2 : [ac*tg] [--*--]
#C1 : [ac*tg]  [--*--]
#C2 : [ac*t+g] [--*-+-]

#C1 : [ac*tg]  [-----]
#C2 : [ac*t+g] [----+-]

sub SIR
{
	my ($t1, $t2) = @_;
	my $transform = Finishing::Assembly::Transform->new;
	my @t1ptu = @{$t1->{padded_to_unpadded}};
	my @t1utp = @{$t1->{unpadded_to_padded}};
	my @t2ptu = @{$t2->{padded_to_unpadded}};
	my @t2utp = @{$t2->{unpadded_to_padded}};
	my @t3ptu = @{$t1->{padded_to_unpadded}};
	my @t3utp = @{$t1->{unpadded_to_padded}};
	
	my $pad_char = $t2->pad_char;
	my $last_num = 0;
	#first we fix the padded array
	for(my $i=0;$i<@t3ptu;$i++)
	{
		if($t3ptu[$i] eq $pad_char)
		{
			if($last_num>0)
			{
				splice(@t3ptu,$i,1);print $last_num."\n";
				splice(@t3ptu,$last_num+1,0,($t3ptu[$last_num]+1));
				
			}
			else
			{
				$t3ptu[$i] = 0;
			}
			$last_num++;
			for(my $j=$last_num+1;$j<@t3ptu;$j++) { $t3ptu[$j] +=1 if($t3ptu[$j] =~ /\d+/); }
		}
		elsif($t3ptu[$i] =~ /\d+/)
		{
			$last_num = $i;		
		}	
	}
	
	#next we fix the unpadded array
	@t3utp = @t3ptu;
	for(my $i=0;$i<@t3utp;$i++) { $t3utp[$i] = $i if($t3utp[$i]=~/\d+/); }
	for(my $i=(scalar(@t3utp)-1);$i>=0;$i--)
	{
		unless($t3utp[$i] =~ /\d+/)
		{
			splice(@t3utp,$i,1);		
		}
	}
	$transform->{unpadded_to_padded} = \@t3utp;
	$transform->{padded_to_unpadded} = \@t3ptu;
	return $transform;
}


sub combine_pads
{
	my ($t1, $t2) = @_;
	my $transform = Finishing::Assembly::Transform->new;
	my @t1ptu = @{$t1->{padded_to_unpadded}};
	my @t1utp = @{$t1->{unpadded_to_padded}};
	my @t2ptu = @{$t2->{padded_to_unpadded}};
	my @t2utp = @{$t2->{unpadded_to_padded}};
	my @t3utp;
	my @t3ptu;
	#the easiest way to do this is to build a list of insertions
	#for each transform
	my @combined_ins_list;
	
	for(my $i = 1;$i<@t2utp;$i++)
	{
		if(($t2utp[$i]-$t2utp[$i-1])>1)
		{
			$combined_ins_list[$i-1] = [$i-1, join('', @t2ptu[ ($t2utp[$i-1]+1)..($t2utp[$i]-1) ]) ];
		}
		if(($t1utp[$i]-$t1utp[$i-1])>1)
		{
			if(!defined $combined_ins_list[$i-1])
			{
				$combined_ins_list[$i-1] =  [$i-1, join('', @t1ptu[ ($t1utp[$i-1]+1)..($t1utp[$i]-1) ]) ];
			}
			else
			{
				$combined_ins_list[$i-1][1] .= join('', @t1ptu[ ($t1utp[$i-1]+1)..($t1utp[$i]-1) ]);
			}
		}
	}
	{
		my @temp;
		foreach(@combined_ins_list) 
		{ 
			push (@temp, $_) if (defined $_);
		} 
		@combined_ins_list = @temp;
	}
	
	@combined_ins_list = reverse sort { cmp_ins($a, $b) } @combined_ins_list;
	my @padded_array;for(my $i=0;$i<@t1utp;$i++){$padded_array[$i] = $i;}

	foreach(@combined_ins_list)
	{
		splice(@padded_array,$_->[0]+1,0,split (//,$_->[1]));

	}
	
	my @unpadded_array = @padded_array;
	for(my $i=0;$i<@unpadded_array;$i++)
	{
		$unpadded_array[$i] = $i if($unpadded_array[$i] =~ /\d+/);
	}
	for(my $i=(scalar (@padded_array)-1);$i>=0;$i--)
	{
		unless( $unpadded_array[$i] =~ /\d+/)
		{
			splice(@unpadded_array,$i,1);
		}
	}
	
	
	$transform->{unpadded_to_padded} = \@unpadded_array;
	$transform->{padded_to_unpadded} = \@padded_array;
	
	return $transform;
}

sub get_pad_list
{
	my ($t1) = @_;
	my @t1ptu = @{$t1->{padded_to_unpadded}};
	my @t1utp = @{$t1->{unpadded_to_padded}};
	
	my %pad_list;
	for(my $i = 1;$i<@t1utp;$i++)
	{
		if(($t1utp[$i]-$t1utp[$i-1])>1)
		{
			my @temp_string = @t1ptu[ ($t1utp[$i-1]+1)..($t1utp[$i]-1) ];
			my %char_hash;
			for(my $j=0;$j<@temp_string;$j++)
			{
				if(exists $char_hash{$temp_string[$j]})
				{
					$char_hash{$temp_string[$j]}++;
				}
				else
				{
					$char_hash{$temp_string[$j]} = 1;
				}				
			}
			$pad_list{$i} = \%char_hash;
		}
	}
	
	return \%pad_list;
}

sub compact_pads
{
	my ($t1) = @_;
	my $transform = Finishing::Assembly::Transform->new;
	my @t1ptu = @{$t1->{padded_to_unpadded}};
	my @t1utp = @{$t1->{unpadded_to_padded}};
	my @t2ptu;
	my @t2utp;
	
	my $pad_list = get_pad_list($t1);
	my @new_pad_list;
	foreach(keys %{$pad_list})#produces a compacted pad list that is the max of the inserted pad values.
	{
		push(@new_pad_list, [$_, max(values %{$pad_list->{$_}})]);
	}
	
	for(my $i=0;$i<@t1utp;$i++){$t2utp[$i] = $t2ptu[$i] = $i;}
	@new_pad_list = reverse sort { cmp_ins($a, $b) } @new_pad_list;
	foreach(@new_pad_list)
	{
		splice(@t2ptu,$_->[0],0, '*'x $_->[1]); #split (//,$_->[1]));
		for(my $i=$_->[0];$i<@t2utp;$i++) { $t2utp[$i] += $_->[1] if($t2utp[$i] =~ /\d+/); }
	}
	
	$transform->{unpadded_to_padded} = \@t2utp;
	$transform->{padded_to_unpadded} = \@t2ptu;
	
	return $transform;		
}

1;
