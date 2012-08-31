package Genome::Model::Tools::Pcap::Sources::Ace::Contig;
our $VERSION = 0.01;

use strict;
use warnings;
use Carp;


use Genome::Model::Tools::Pcap::Sources::Ace::Read;
use Storable;
use base (qw(Genome::Model::Tools::Pcap::Sources::Ace::SequenceItem));

=cut
Contig Data Structure
Contig:
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
    
    return $self;   
}

sub get_map {
    my ($self) = @_;
    return { self => [ 'self' ],
			 contig   => [ 'name', 'complemented', 'base_count', 'length', 'read_count', 'base_seg_count' ],
             quality => [ 'unpadded_base_quality' ],
             base_segments => [ 'base_segments' ],
             tags => [ 'tags' ],
             children => [ 'children' ],
   			 padded_base_string   => [ 'padded_base_string', 'unpadded_base_string', 'alignment' ],
			 unpadded_base_string => [ 'unpadded_base_string', 'padded_base_string' ],
			 padded_base_quality   => [ 'padded_base_quality', 'unpadded_base_quality', 'alignment' ],
			 unpadded_base_quality => [ 'unpadded_base_quality', 'padded_base_quality' ]
           };
}

sub padded_base_string
{
    my ($self,$object) = @_;
    return 1 unless (@_ > 1);
    
    my $consensus;
    my $fh = $self->{fh};
    my $index = $self->{index};
    $fh->seek($index->{offset},0);
    
    my $line = <$fh>;	
	#$fh->read($consensus, 1000);
    while ($line = <$fh>) 
    {
        if (length($line) == 1) 
        {
            last;
        }
		if ( (substr( $_, 0, 1 ) eq " ") ||(substr( $_, 0, 1 ) eq "\t" )) 
		{
        	last if ( $_ =~ /^\s*$/ );
        }		
		
        chomp $line;        
        $consensus .= $line;
    }
	$object->{just_load} = 1;
    $object->padded_base_string($consensus);
	$object->{just_load} = 0;
    return 1; 
}

sub padded_base_quality
{
    return 0;           
}

sub unpadded_base_string
{
    return 0;    
}

sub unpadded_base_quality
{
    my ($self,$object) = @_;
    return 1 unless (@_ > 1);
    my $fh = $self->{fh};
    my $index = $self->{index};
    $fh->seek($index->{base_qualities}{offset},0);       
    while (my $line = <$fh>) 
    {
        if ($line =~ /^BQ/) 
        {
            last;
        }
    }
    my @bq;
    while (my $line = <$fh>) 
    {
        if ($line =~ /^\s*$/) 
        {
            last;
        }
        chomp $line;
        $line =~ s/^ //; # get rid of leading space
        push @bq, split(/ /,$line);
    }
	$object->{just_load} = 1;
    $object->unpadded_base_quality(\@bq);  
	$object->{just_load} = 0;
    return 1; 
}

sub has_alignment
{
    return "padded_base_string";
}

#methods inherited from Item Source

sub position 
{
    my ($self,$object) = @_;    
    return 0;    
}


sub length 
{
    my ($self,$object) = @_;
    return 1 unless (@_ > 1);
    my $fh = $self->{fh};
    my $index = $self->_index;
    my $reader = $self->{reader};
    my $name = $self->{name};
    $fh->seek($index->{offset},0);
    my $line = <$fh>;
    chomp $line;
    my @tokens = <$fh>;
    my %ret_val = (
        base_count => $tokens[2],
        read_count => $tokens[3],
        base_seg_count => $tokens[4],        
        u_or_c => $tokens[5],
    );
	$object->{just_load} = 1;    
    $object->length($ret_val{base_count});
    $object->base_count($ret_val{base_count}) if !$object->already_loaded("base_count");
    $object->read_count($ret_val{read_count}) if !$object->already_loaded("read_count");
    $object->base_segment_count($ret_val{base_seg_count}) if !$object->already_loaded("base_segment_count");
    $object->complemented(lc($ret_val{u_or_c})eq"C"?1:0) if !$object->already_loaded("complemented");
	$object->{just_load} = 0;
    return 1;   
}

sub _build_read_tag {
    my ($self, $obj) = @_;
    my $tag = new Genome::Model::Tools::Pcap::Tag(
        type => $obj->{tag_type},
        date => $obj->{date},
        source => $obj->{program},
        parent => $obj->{read_name},
        scope => 'ACE',
        start => $obj->{start_pos},
        stop => $obj->{end_pos},
    );
    return $tag;
}


sub tags
{
    my ($self,$object) = @_;
    return 1 unless (@_ > 1);  
    my $input = $self->{fh};
    my $index = $self->_index;
    my $reader = $self->{reader};
    my $name = $self->{name};
    my @tags = @{$index->{contig_tags}};  
    my @contig_tags;
    foreach my $tag_index (@tags)
    {
        $input->seek($tag_index->{offset},0);
        my $contig_tag = Genome::Model::Tools::Pcap::TagParser->new()->parse($input);
        push @contig_tags, $contig_tag;    
    }
	$object->{just_load} = 1;
    $object->tags( \@contig_tags );
	$object->{just_load} = 0;
    return 1;
}

sub copy
{
    my ($self) = @_;
    return dclone ($self);
}

sub copy_tag 
{
    my ($self, $tag) = @_;
    return dclone $tag;
}

sub start_position
{
    my ($self,$object) = @_;
    return 1 unless (@_ > 1);
	$object->{just_load} = 1;
    $object->start_position($self->position);
	$object->{just_load} = 0;
    return 1;    
}

sub end_position 
{
    my ($self,$object) = @_;
    return 1 unless (@_ > 1);
	$object->{just_load} = 1;
    $object->end_position($self->position + $self->length);
	$object->{just_load} = 0;
    return 1;
}

sub base_count
{
    my ($self,$object) = @_;
    return 1 unless (@_ > 1);
    my $fh = $self->{fh};
    my $index = $self->_index;
    my $reader = $self->{reader};
    my $name = $self->{name};
    $fh->seek($index->{offset},0);
    my $line = <$fh>;
    chomp $line;
    my @tokens = split / /, $line;
    my %ret_val = (        
        base_count => $tokens[2],
        read_count => $tokens[3],
        base_seg_count => $tokens[4],        
        u_or_c => $tokens[5],
    );
	$object->{just_load} = 1;         
    $object->base_count($ret_val{base_count});#probably should add length
    $object->read_count($ret_val{read_count}) if !$object->already_loaded("read_count"); 
    #$object->base_segment_count($ret_val{base_seg_count}) if !$object->already_loaded("base_segment_count");
    $object->complemented(lc($ret_val{u_or_c})eq"C"?1:0) if !$object->already_loaded("complemented");
	$object->{just_load} = 0;
    return 1;       
}

sub base_segment_count
{
    my ($self,$object) = @_;
    return 1 unless (@_ > 1);
	$object->{just_load} = 1;
    $object->base_segment_count($self->_index->{base_segments}{line_count});
	$object->{just_load} = 0;    
    return 1;
}

sub complemented
{
    my ($self,$object) = @_;
    return 1 unless (@_ > 1);
    my $fh = $self->{fh};
    my $index = $self->_index;
    my $reader = $self->{reader};
    my $name = $self->{name};
    $fh->seek($index->{offset},0);
    my $line = <$fh>;
    chomp $line;
    my @tokens = split / /,$line;
    my %ret_val = (        
        base_count => $tokens[2],
        read_count => $tokens[3],
        base_seg_count => $tokens[4],        
        u_or_c => $tokens[5],
    );
	$object->{just_load} = 1;    
    $object->base_count($ret_val{base_count}) if !$object->already_loaded("base_count");
    $object->read_count($ret_val{read_count}) if !$object->already_loaded("read_count"); 
    #$object->base_segment_count($ret_val{base_seg_count}) if !$object->already_loaded("base_segment_count");
    $object->complemented(lc($ret_val{u_or_c})eq"C"?1:0);
	$object->{just_load} = 0;
    return 1;      
    
}

sub read_count
{
    my ($self,$object) = @_;
    return 1 unless (@_ > 1);
    my $fh = $self->{fh};
    my $index = $self->_index;
    my $reader = $self->{reader};
    my $name = $self->{name};
    $fh->seek($index->{offset},0);
    my $line = <$fh>;
    chomp $line;
    my @tokens = split / /,$line;
    my %ret_val = (        
        base_count => $tokens[2],
        read_count => $tokens[3],
        base_seg_count => $tokens[4],        
        u_or_c => $tokens[5],
    );
	$object->{just_load} = 1;    
    $object->base_count($ret_val{base_count}) if !$object->already_loaded("base_count");
    $object->read_count($ret_val{read_count});
    #$object->base_segment_count($ret_val{base_seg_count}) if !$object->already_loaded("base_segment_count"); 
    $object->complemented(lc($ret_val{u_or_c})eq"C"?1:0) if !$object->already_loaded("complemented");
	$object->{just_load} = 0;
    return 1;
}
    

sub base_segments
{
    my ($self,$object) = @_;
    return 1 unless (@_ > 1);
    my $fh = $self->{fh};
    my $index = $self->_index;
    my $reader = $self->{reader};
    $fh->seek($index->{base_segments}{offset},0);
    my @base_segments;
    while(my $obj = $reader->next_object)
    {
        last if($obj->{type} ne "base_segment");
        push @base_segments, $obj;
    }
	$object->{just_load} = 1;
    $object->base_segments(\@base_segments);    
	$object->{just_load} = 0;
    return 1;
}

sub children
{
    my ($self,$object) = @_;
    return 1 unless (@_ > 1);
    my %children;
    foreach my $read_index (values %{$self->_index->{reads}} )
    {
        my $read_callback = Genome::Model::Tools::Pcap::Sources::Ace::Read->new(name => $read_index->{read}{name},
        index => $self->_index, reader => $self->{reader} );
        my $read = Genome::Model::Tools::Pcap::Read->new(callbacks => $read_callback);
        $children{$read_index->{read}{name}} = $read;        
    }
	$object->{just_load} = 1;
    $object->children(\%children);
	$object->{just_load} = 0;
    return 1;
}

sub _index
{
	my ($self) = @_;
    return $self->{index};

}




1;


