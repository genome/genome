#The Contig Merging and Splitting Utility Library has several classes.

package Genome::Model::Tools::Pcap::Transform;
our $VERSION = 0.01;

use strict;
use warnings;
use Carp;
   
my $pkg = "Genome::Model::Tools::Pcap::Transform";
sub new#($padded_string, $pad_char)#return $transform
{
	croak("$pkg:new:no class given, quitting") if @_ < 1;
    
	my ($caller, $padded_string, $pad_char) = @_;
	my $caller_is_obj = ref($caller);
    my $class = $caller_is_obj || $caller;
    $pad_char = '*' if (@_<3);
	my $transform = {};
	bless ($transform, $class);
	$transform->{offset}=0;
    
    if(@_ == 1)
    {
        return $transform;   
    }
	
	my $temp_length = length($padded_string);
	my @padded_to_unpadded;
	my @unpadded_to_padded;
	my $unpadded_index = 0;
    $#padded_to_unpadded = $temp_length-1;
	for(my $i = 0;$i < $temp_length;$i++)
	{
		$unpadded_to_padded[$unpadded_index]=$i;
        if(substr($padded_string, $i, 1) ne $pad_char)
		{
			$padded_to_unpadded[$i]=$unpadded_index;
			$unpadded_index++;
		}
		else
		{
			$padded_to_unpadded[$i]='*';
		}

	}
	$transform->{padded_to_unpadded} = \@padded_to_unpadded;
	$transform->{unpadded_to_padded} = \@unpadded_to_padded;
	$transform->{pad_char} = $pad_char;
	return $transform;
}

sub loaded
{
    my ($transform) = @_;
	if(defined $transform->{padded_to_unpadded})
	{
		return 1;
	}
	else
	{
		return 0;
	}
}

sub _offset
{
	my ($self,$offset) = @_;
	$self->{offset} = $offset if @_>1;
	return $self->{offset};
}

sub derive_from_base_string
{
	my ($self, $padded_string,$pad_char) = @_;
	$pad_char = '*' if (@_<3);
	my $temp_length = length($padded_string);
	my @padded_to_unpadded;
	my @unpadded_to_padded;
	my $unpadded_index = 0;
    $#padded_to_unpadded = $temp_length-1;
	for(my $i = 0;$i < $temp_length;$i++)
	{
		$unpadded_to_padded[$unpadded_index]=$i;
        if(substr($padded_string, $i, 1) ne $pad_char)
		{
			$padded_to_unpadded[$i]=$unpadded_index;
			$unpadded_index++;
		}
		else
		{
			$padded_to_unpadded[$i]='*';
		}

	}
	$self->{padded_to_unpadded} = \@padded_to_unpadded;
	$self->{unpadded_to_padded} = \@unpadded_to_padded;
	$self->{pad_char} = $pad_char;	

}
sub derive_from_quality_array
{
	my ($self, $padded_array,$pad_char) = @_;
	$pad_char = '*' if (@_<3);
	my $temp_length = scalar @{$padded_array};
	my @padded_to_unpadded;
	my @unpadded_to_padded;
	my $unpadded_index = 0;
    $#padded_to_unpadded = $temp_length-1;
	for(my $i = 0;$i < $temp_length;$i++)
	{
		$unpadded_to_padded[$unpadded_index]=$i;
        if($$padded_array[$i] ne $pad_char)
		{
			$padded_to_unpadded[$i]=$unpadded_index;
			$unpadded_index++;
		}
		else
		{
			$padded_to_unpadded[$i]='*';
		}
	}
	$self->{padded_to_unpadded} = \@padded_to_unpadded;
	$self->{unpadded_to_padded} = \@unpadded_to_padded;
	$self->{pad_char} = $pad_char;
	return $self;

}
sub pad_string#($string, $pad_char)#returns $out_string
{
	my @out_string;	
	my ($transform, $string, $pad_char) = @_;
	my $use_pad_char = 0;
	if(@_ > 2)
	{
		$use_pad_char = 1;
	}
	if(@{$transform->{unpadded_to_padded}} != length ($string))
	{
		print "padded_to_unpadded is not equal to length of $string\n";
	}
	if($use_pad_char)
	{
		for(my $i=0;$i<@{$transform->{padded_to_unpadded}};$i++)
		{
			$out_string[$i] = $pad_char;#this will preserve pad chars, in case we
			#are doing a composite transform
		}
	}
	else
	{
		for(my $i=0;$i< @{$transform->{padded_to_unpadded}};$i++)
		{
			$out_string[$i] = $transform->{padded_to_unpadded}[$i];#this will preserve pad chars, in case we
			#are doing a composite transform
		}
	}
	my $temp_length = length($string);
	for(my $i=0;$i<$temp_length;$i++)
	{
		$out_string[$transform->{unpadded_to_padded}[$i]]= substr($string, $i, 1);
	}

	my $out_string = join '', @out_string;	
	return $out_string;
}
#this function will pad a string that only overlaps the region to be padded, specify an offset
#and it will pad up to the length of string
sub pad_string_partial#($string, $offset, $pad_char)#returns $out_string
{
	my @out_string;	
	my ($transform, $string, $offset, $pad_char) = @_;
	my $use_pad_char = 0;
	my $pre_string; #string that comes before the merge sequence
	my $post_string; #string that over_hangs the end
	if(@_ > 3)
	{
		$use_pad_char = 1;
	}

	if($use_pad_char)
	{
		for(my $i=0;$i<@{$transform->{padded_to_unpadded}};$i++)
		{
			$out_string[$i] = $pad_char;
		}
	}
	else
	{
		for(my $i=0;$i<@{$transform->{padded_to_unpadded}};$i++)
		{
			$out_string[$i] = $transform->{padded_to_unpadded}[$i];#this will preserve pad chars, in case we
			#are doing a composite transform
		}
	}
	my $temp_length = length($string);
	#fix need to account for strings that start before or after consensus
	my $over_hang = 0;#this is the amount that hangs off the end of the contig
	my $pre_hang = 0;
	if(($temp_length + $offset)>@{$transform->{unpadded_to_padded}})
	{
		$over_hang = $temp_length;
		$temp_length = @{$transform->{unpadded_to_padded}}-$offset;
		$over_hang -= $temp_length;
		$post_string = substr($string, $temp_length,$over_hang); 	
	}
	if($offset<0)
	{
		$pre_hang = $offset;
		#$offset=0;
		$pre_string = substr( $string, 0, -$pre_hang);
		#$string = substr($string, -$pre_hang,length($string)+$pre_hang);
	}


	for(my $i=0;$i<$temp_length;$i++)
	{
		next if (($i+$offset)<0);
		$out_string[$transform->{unpadded_to_padded}[$i+$offset]]= substr($string, $i, 1);
	}
	#clip off first and last parts of outstring if necessary
	#the reason we add a 1 at the end is because we are getting the length of the padded sequence inclusuve of 
	#the last base, NOT the last base+pads, i.e. we want ACT**G, NOT ACT**G**, which is what we would get if
	#added the 1 to temp_length, i.e. 1 subtracted from ACT**G**T is ACT**G**,
	if($offset<0)
	{
		@out_string = splice(@out_string, $transform->{unpadded_to_padded}[0], $transform->{unpadded_to_padded}[$temp_length+$offset-1]-$transform->{unpadded_to_padded}[0]+1);
	}
	else
	{
		@out_string = splice(@out_string, $transform->{unpadded_to_padded}[$offset], $transform->{unpadded_to_padded}[$temp_length+$offset-1]-$transform->{unpadded_to_padded}[$offset]+1);
	}
	my $out_string = join '', @out_string;
	if($pre_hang<0)
	{
		$out_string = $pre_string.$out_string;
	}
	if($over_hang>0)
	{
		$out_string = $out_string.$post_string;
	}	
	return $out_string;
}
sub unpad_string#($string, $pad_char )#returns $out_string
{
	my ($transform, $string, $pad_char) = @_;
	if(@_ < 3)
	{
		$pad_char = $transform->{pad_char};
	}
	my $temp_length = length ($string);
	my $out_string;
	for(my $i=0;$i<$temp_length;$i++)
	{
		if($transform->{padded_to_unpadded}[$i]=~/\d+/)#if it's a numerical position
		{
			$out_string .= substr($string, $i, 1);
		}
	}
	return $out_string;

}

# given the following three strings
# 1.           SSSSSSSSSSSSSSAS
# 2.  SSSSSSSSSSSSSSSSSSSSSSSASSSSSSSSSS
# 3.  SSSSSSSSSSSSSSSSSSSSS*SSASSSS*SSSSSS
# Where 1 is a string that starts at $offset from the beginning of string2, and string 2 is an unpadded
# version of string 3, this function takes an $transform, $nUnpadPosition (a position in string 1), and
# $offset, which is the offset of string 1 from the start of string 2, and returns the coordinate of this
# position in string 3 as output.  i.e., it returns the position of A in string 3, given the position of A
# in string 1, and the offset of the beginning of string 1 to the beginning of string 2.  
# This function is convenient for taking a read index, and finding it's position in the new consensus that is
# produced by bioperl.  
sub get_pad_position#($unpadded_position, $offset)#returns Padded Position in string 3
{
	my @out_string;	
	my ($transform, $unpadded_position, $offset) = @_;
	
	if(@_<3)
	{
		$offset=0;
	}
	if($unpadded_position < 1)
	{
		return $unpadded_position+$transform->_offset;	
	}
	elsif($unpadded_position > scalar @{$transform->{unpadded_to_padded}})
	{
		my $tempoffset = $unpadded_position - scalar @{$transform->{unpadded_to_padded}};
		return @{$transform->{padded_to_unpadded}}+$tempoffset+$transform->_offset;	
	}
	my $pad_string_position = $unpadded_position + $offset;
	if($pad_string_position<1) 
	{
		return $pad_string_position;
	}
	elsif($pad_string_position>@{$transform->{unpadded_to_padded}})
	{
		return $pad_string_position-@{$transform->{unpadded_to_padded}}+@{$transform->{padded_to_unpadded}};
	}
	return ${$transform->{unpadded_to_padded}}[$pad_string_position-1]+1+$transform->_offset;		
}

sub get_unpad_position
{
	my ($transform, $padded_position) = @_;
	if($padded_position < 1)
	{
		return $padded_position+$transform->_offset;	
	}
	elsif($padded_position > scalar @{$transform->{padded_to_unpadded}})
	{
		my $tempoffset = $padded_position - scalar @{$transform->{padded_to_unpadded}};
		return @{$transform->{unpadded_to_padded}}+$tempoffset+$transform->_offset;		
	}
	while($padded_position>1)
	{
		if(${$transform->{padded_to_unpadded}}[$padded_position-1] eq $transform->{pad_char})
		{
			$padded_position--;
		}
		else
		{
			last;
		}
	}
	#if we've hit 0 and it's still a pad, return 0
	return 1 if($padded_position==1 and ${$transform->{padded_to_unpadded}}[$padded_position-1] eq $transform->{pad_char});
	
	return ${$transform->{padded_to_unpadded}}[$padded_position-1]+1+$transform->_offset;
}

sub pad_array
{
    my ($self, $unpadded_array, $pad_char) = @_;
    
    my @out_array; 
    #my ($transform, $string, $pad_char) = @_;
    my $use_pad_char = 0;
    if(@_ > 2)
    {
        $use_pad_char = 1;
    }
    if(@{$self->{unpadded_to_padded}} != @{$unpadded_array})
    {
        print "padded_to_unpadded is not equal to length of array \n";
    }
    if($use_pad_char)
    {
        for(my $i=0;$i<@{$self->{padded_to_unpadded}};$i++)
        {
            $out_array[$i] = $pad_char;#this will preserve pad chars, in case we
            #are doing a composite transform
        }
    }
    else
    {
        for(my $i=0;$i< @{$self->{padded_to_unpadded}};$i++)
        {
            $out_array[$i] = $self->{padded_to_unpadded}[$i];#this will preserve pad chars, in case we
            #are doing a composite transform
        }
    }
    my $temp_length = @{$unpadded_array};;
    for(my $i=0;$i<$temp_length;$i++)
    {
        $out_array[$self->{unpadded_to_padded}[$i]]= $$unpadded_array[$i];
    }

    return \ @out_array;
}

sub unpad_array
{
    my ($self, $padded_array, $pad_char) = @_;
    if(@_ < 3)
    {
        $pad_char = $self->{pad_char};
    }
    my $temp_length = @{$padded_array};
    my @out_array;
    
    for(my $i=0;$i<$temp_length;$i++)
    {
        if($self->{padded_to_unpadded}[$i]=~/\d+/)#if it's a numerical position
        {
            push @out_array, $$padded_array[$i];
        }
    }
    return \ @out_array;
}

sub unpad_length
{
	my ($self) = @_;
	return scalar @{$self->{unpadded_to_padded}};
}

sub pad_length
{
	my ($self) = @_;
	return scalar @{$self->{padded_to_unpadded}};
}


sub copy
{
    my ($self) = @_;    
    croak("$pkg:new:no class given, quitting") if @_ < 1;
    my ($caller, $padded_string, $pad_char) = @_;
    my $caller_is_obj = ref($caller);
    my $class = $caller_is_obj || $caller;    
    my $copy = {unpadded_to_padded => [ @{$self->{unpadded_to_padded}}],
                 padded_to_unpadded => [@{$self->{padded_to_unpadded}}],
                 pad_char => $self->{pad_char},
				 offset => $self->{offset}};
    bless ($copy, $class);
    
    return $copy;
}
1;

