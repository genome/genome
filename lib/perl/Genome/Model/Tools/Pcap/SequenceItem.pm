package Genome::Model::Tools::Pcap::SequenceItem;
our $VERSION = 0.01;

use strict;
use warnings;
use Carp;
use Genome::Model::Tools::Pcap::Tag;
use Genome::Model::Tools::Pcap::Item;
#use Genome::Model::Tools::Pcap::Transform;
use base (qw(Genome::Model::Tools::Pcap::Item));

sub new {
    croak("__PACKAGE__:new:no class given, quitting") if @_ < 1;
    my ($caller, %params) = @_;
    my $caller_is_obj = ref($caller);
    my $class = $caller_is_obj || $caller;
    my $self = $class->SUPER::new(%params);
    #my $self = {};# = {%arg};
    bless ($self, $class);    
    $self->{recent_bases} = 'both';
	$self->{recent_quals} = 'both';
	$self->{recent_chrom} = 'both';
	$self->{recent_align} = 'both';
	$self->{always_update} = 1;
	if(exists $params{callbacks})
	{
		$self->callbacks($params{callbacks});
	}
    return $self;
}

sub get_base_at_position
{
    my ($self, $position) = @_;

    $self->fatal_msg("You must supply a base position to get that base") unless
	$position =~ /^\d+$/;

    my $bases = $self->padded_base_string;
    
    return substr($bases, $position,1);
}

sub length #no setter for length, since it doesn't make sense to allow someone to set the length for an item that
#already has a length (i.e.  if you have a string ACTG, it doesn't make sense to set the length to anything other
#than 4, which is already implied
{
    my ($self,$type) =@_;
	
	if($self->_has_alignment)
    {
		return length($self->unpadded_base_string) if (defined $type && $type eq "unpadded");
        return length($self->padded_base_string);
    }
    elsif($self->unpadded_base_string)
    {
		warn "Sequence does not have alignment information!\n" if (defined $type && $type eq "padded");
        return length($self->unpadded_base_string);
    }
    else #invalid
    {
        return 0;
    }
}

sub length_deprecated
{
	my ($self, $type) = @_;
	
	if(defined $type)
	{
		if($type eq "unpadded")
		{
			return $self->_transform->unpad_length;									
		}
		else
		{
			return $self->_transform->pad_length;
		}	
	}
}

sub freeze
{
	my ($self) = @_;
	$self->SUPER::freeze;
	#if($self->already_loaded("sequence"))
	#{
	#	$self->freeze;
	#}
}

sub thaw
{
	my ($self)  = @_;
	$self->SUPER::thaw;
	#if($self->already_loaded("sequence"))
	#{
	#	$self->thaw($obj,$file_name, $fh);
	#}
}

# sequence methods
sub _transform
{
    my ($self, $transform) = @_;
    
    $self->{transform} = $transform if (@_ > 1);
    
    if($self->{recent_align} eq "transform")
	{
		return $self->{transform};
	}
	else
	{
		$self->{transform} = Genome::Model::Tools::Pcap::Transform->new if(!defined $self->{transform});
		$self->_load_transform;
		$self->{recent_align} = "transform";
        return $self->{transform};
	}
}

sub get_transform
{
    my ($self) = @_;
	$self->load_data("padded_base_string");
	return $self->_transform->copy($self->_transform);    
}

sub _load_transform
{
    my ($self) = @_;
	
	if(my $string = $self->check_and_load_data("padded_base_string"))
	{
		$self->{transform}->derive_from_base_string($string,'*');
	}		
	else
	{
		die "Could not derive Alignment Information!\n";
	}	
}

sub _has_alignment
{
	my ($self) = @_;
	if($self->{transform}||$self->padded_base_string)
	{
		return 1;
	}
	else
	{
		return 0;
	}
}

sub padded_base_string
{
    my ($self, $value) = @_;
    my ($name) = (caller(0))[3] =~ /.+::(.+)/;
    if(@_ >1)
    {
		#this is a hack, need to add a mechanism for verification of transforms after they have change,
		#until then, the only solution for updating transforms is to nullify the padded versions of
		#the quality arrays
		$self->{recent_quals} = 'unpadded_base_quality' if($self->{recent_quals} =~/both/);
		if($self->{recent_quals}=~/^padded_base_quality/)#should probably check recent_quals instead
		{$self->check_and_load_data("unpadded_base_quality",$self->check_and_load_data("unpadded_base_quality"));}
#		if(defined $self->{padded_chromat_positions})
#		{$self->check_and_load_data("unpadded_chromat_positions",$self->check_and_load_data("unpadded_chromat_positions"));}		
		$self->check_and_load_data($name, $value);
        $self->{recent_bases} = $name;
		$self->{recent_align} = $name;
		
				
    }  
	my $string;  
    if(($self->{recent_bases} =~ /^$name|both/) && ($string = $self->check_and_load_data($name)))
	{
		return $string;
	}
	elsif(($self->{recent_bases} =~ /^$name|both/)&& 
		  ($string = $self->check_and_load_data("unpadded_base_string"))&&
		  $self->_has_alignment) 
	{
		$self->{recent_bases} = 'both';
		$self->{just_load} = 1;#don't want to register a derivation as a data change
		my $temp = $self->check_and_load_data("padded_base_string", $self->_transform->pad_string($string));
		$self->{just_load} = 0;
		return $temp;
    	#return $self->_transform->pad_string($string);		
    }
	elsif(!$self->{always_update}&& ($string = $self->check_and_load_data($name)))
	{
		return $string;
	}
    else
    {
 		warn "Data is undefined\n";
        return "";
    }    
}

sub padded_base_quality
{
    my ($self, $value) = @_;
    my ($name) = (caller(0))[3] =~ /.+::(.+)/;
    if(@_>1)
    {        	    
        $self->check_and_load_data($name, $value);
		$self->check_and_load_data("unpadded_base_quality", []);
		$self->{recent_quals} = $name;                           
    }
	my $base_qual;
	if($self->{recent_quals} =~ /^$name|both/ && 
	   ($base_qual = $self->check_and_load_data($name)))
	{
		return $base_qual;
	}
	elsif($self->{recent_quals} =~ /unpadded_base_quality|both/ && 
		  ($base_qual = $self->check_and_load_data("unpadded_base_quality")) )
	{
		$self->{recent_quals} = 'both';
		$self->{just_load} = 1;#don't want to register a derivation as a data change
		my $temp = $self->check_and_load_data("padded_base_quality", $self->_transform->pad_array($base_qual));
    	$self->{just_load} = 0;
		return $temp;
		#return $self->_transform->pad_array($base_qual);		
    }
	elsif(!$self->{always_update}&&
		  ($base_qual = $self->check_and_load_data($name)))
	{
		return $base_qual;
	}
    else
    {
		die "Data is undefined\n";
        return [];
    }        
}

sub padded_chromat_positions
{
	my ($self, $value) = @_;
    my ($name) = (caller(0))[3] =~ /.+::(.+)/;
    if(@_>1)
    {        	    
        $self->check_and_load_data($name, $value);
		$self->check_and_load_data("unpadded_chromat_positions", undef);
		$self->{recent_chrom} = $name;                           
    }
	my $chrom;
	if($self->{recent_chrom} =~ /^$name|both/ && 
	   ($chrom = $self->check_and_load_data($name)))
	{
		return $chrom;
	}
	elsif($self->{recent_chrom} =~ /unpadded_chromat_positions|both/ && 
	     ($chrom = $self->check_and_load_data("unpadded_chromat_positions"))&&
		 $self->_has_alignment )
	{
		$self->{recent_chrom} = 'both';
    	return $self->_transform->pad_array($chrom);
    }
	elsif($self->{always_update} && ($chrom = $self->check_and_load_data($name)))
	{
		return $chrom;
	}
    else
    {
		die "$name is undefined\n";
        return [];
    }	
}


sub unpadded_base_string
{
    my ($self, $value) = @_;
    my ($name) = (caller(0))[3] =~ /.+::(.+)/;
    if(@_ > 1)
    {
        $self->check_and_load_data($name,$value);
		$self->_transform if(defined $self->check_and_load_data("padded_base_string")); #go ahead and derive and save transform info if it's there
		$self->check_and_load_data("padded_base_string",undef);
		$self->{recent_bases} = $name;
    }
	my $string;
    if($self->{recent_bases} =~ /^$name|both/ && ($string = $self->check_and_load_data("unpadded_base_string")))
    {
        return $string;        
    }
    elsif($self->{recent_bases} =~ /padded_base_string|both/ && ($string = $self->check_and_load_data("padded_base_string")))
    {
		$self->{recent_bases} = 'both';
		$self->{just_load} = 1;#don't want to register a derivation as a data change
		my $temp = return $self->check_and_load_data("unpadded_base_string", $self->_transform->unpad_string($string));
        $self->{just_load} = 0;
		return $temp;
		#return $self->_transform->unpad_string($string);        
    }
	elsif(!$self->{always_update}&&($string = $self->check_and_load_data("unpadded_base_string")) )
	{
		return $string;	
	}
    else
    {
		die "$name is undefined\n";
        return '';
    }
}

sub unpadded_base_quality
{
    my ($self, $value) = @_;
    my ($name) = (caller(0))[3] =~ /.+::(.+)/;   
    if(@_ > 1)
    { 
	    $self->check_and_load_data($name,$value);
		#$self->_transform if(defined $self->check_and_load_data("padded_base_quality"));
		#$self->check_and_load_data("padded_base_quality", undef);
		$self->{recent_quals} = $name;
    }
	my $base_qual;
	if($self->{recent_quals} =~ /^$name|both/ &&
	   ($base_qual = $self->check_and_load_data("unpadded_base_quality")))
    {
        return $base_qual;        
    }
    elsif($self->{recent_quals} =~ /padded_base_quality|both/ &&
	      ($base_qual = $self->check_and_load_data("padded_base_quality"))&&
		  $self->_has_alignment)
    {
		$self->{recent_quals} = 'both';
		$self->{just_load} = 1;#don't want to register a derivation as a data change,
		my $temp = $self->check_and_load_data("unpadded_base_quality", $self->_transform->unpad_array($base_qual));
        $self->{just_load} = 0;
		return $temp;
		#return $self->_transform->unpad_array($base_qual);        
    }
	elsif(!$self->{always_update}&&($base_qual = $self->check_and_load_data("unpadded_base_quality")))
	{	
		if($base_qual = $self->check_and_load_data("unpadded_base_quality"))
    	{
    	    return $base_qual;        
    	}		
	}
    else
    {
		die "$name is undefined!\n";
        return [];
    }    
}

sub unpadded_chromat_positions
{
    my ($self, $value) = @_;
    my ($name) = (caller(0))[3] =~ /.+::(.+)/;   
    if(@_ > 1)
    { 
	    $self->check_and_load_data($name,$value);
		#$self->_transform if(defined $self->check_and_load_data("padded_chromat_positions"));
		#$self->check_and_load_data("padded_chromat_positions", undef);
		$self->{recent_chrom} = $name;
    }
	my $chrom;
	if($self->{recent_chrom} =~ /^$name|both/ &&
	   ($chrom = $self->check_and_load_data($name)))
    {
        return $chrom;        
    }
    elsif($self->{recent_chrom} =~ /padded_chromat_positions|both/ &&
	      ($chrom = $self->check_and_load_data("padded_chromat_positions"))&&
		  $self->_has_alignment)
    {
		$self->{recent_chrom} = "both";
        return $self->_transform->unpad_array($chrom);        
    }
	elsif(!$self->{always_update}&&($chrom = $self->check_and_load_data("unpadded_chromat_positions")))
	{	
		if($chrom = $self->check_and_load_data("unpadded_chromat_positions"))
    	{
    	    return $chrom;        
    	}		
	}
    else
    {
		die "$name is undefined!\n";
        return [];
    }    
}

sub get_padded_base_quality
{
	my ($self, $padded_base_position) = @_;
	my $position = ${$self->_transform->{padded_to_unpadded}}[$padded_base_position];
	return 0 if $position eq '*';
	return ${$self->{unpadded_base_quality}}[$position];
}

sub get_padded_base_value
{
	my ($self, $padded_base_position) = @_;
	my $position = ${$self->_transform->{padded_to_unpadded}}[$padded_base_position];
	return "*" if $position eq '*';
	return substr($self->{unpadded_base_string},$position, 1);
}

#new functions

sub replace_chromat_pads
{
	my ($self, $padded_chromat_positions) = @_;
	$padded_chromat_positions = [@{$padded_chromat_positions}];
	my $on_pad = 0;
	my $pad_start_pos = 0;
	my $pad_end_pos = 0;
	for(my $i=0;$i<@{$padded_chromat_positions};$i++)
	{
		if($padded_chromat_positions->[$i] eq '*')
		{
			if($on_pad)
			{
				$pad_end_pos = $i;
			}
			else
			{
				$pad_start_pos = $i;
				$pad_end_pos = $i;
			}
			$on_pad = 1;		
		}
		else 
		{
			if($on_pad)
			{
				$pad_start_pos--;
				$pad_end_pos++;
				my $start_chromat_value = 0;
				my $increment = 0;
				if($pad_start_pos < 0)
				{
					$increment = 5;					
					$start_chromat_value = 5;
					
					
				}
				else
				{
					$increment = ($padded_chromat_positions->[$pad_end_pos] - 
					$padded_chromat_positions->[$pad_start_pos])/($pad_end_pos - $pad_start_pos);
					$start_chromat_value = $padded_chromat_positions->[$pad_start_pos] + $increment;
				}
				for(my $j=($pad_start_pos+1);$j<$pad_end_pos;$j++)
				{
					$padded_chromat_positions->[$j]=int($start_chromat_value + $increment*($j-($pad_start_pos+1)));
				}				
				$on_pad = 0;
			}		
		}	
	}
	if($on_pad)
	{
		$pad_start_pos--;
		$pad_end_pos++;		
		my $start_chromat_value = 0;
		my $increment = 0;
		if($pad_start_pos < 0)
		{				
			$increment = 5;
			$start_chromat_value = 0;
		}
		else
		{
			$increment = 5;
			$start_chromat_value = $padded_chromat_positions->[$pad_start_pos] + $increment;
		}
		for(my $j=($pad_start_pos+1);$j<$pad_end_pos;$j++)
		{
			$padded_chromat_positions->[$j]=int($start_chromat_value + $increment*($j-($pad_start_pos+1)));
		}				
		$on_pad = 0;
	}
	return $padded_chromat_positions;
}

sub padded_base_string_xn_positions
{
    my $self = shift;

    my @tmp = split ('', $self->padded_base_string);
    my @xn_pos;
    my $count = 0;

    foreach my $base (@tmp)
    {
	    $count++;
	    if ($base =~ /^[xn]$/i)
	    {
	        push @xn_pos, $count;
	    }
    }
    return \@xn_pos;
}

sub copy
{
    my ($self,$item) = @_;
    
	return Storable::dclone($item);    
}

1;
