package Genome::Model::Tools::Pcap::Item;
our $VERSION = 0.01;

use strict;
use warnings;
use Carp;
use Storable;
use Genome::Model::Tools::Pcap::Tag;
use Storable;
use base(qw(Genome::Model::Tools::Pcap::Proxy));
my $pkg = "Genome::Model::Tools::Pcap::Item";

=pod

=head1 NAME

Item - Base class for Reads, Contigs, and SuperContigs (through SequenceItem).  This class is not used directly, but instead encapsulates low level functionality that is common to Reads, Contigs, and SuperContigs.

=head1 DESCRIPTION

Genome::Model::Tools::Pcap::Item has a position, length, tags, and children.  The data items can be set through accessor methods.

=head1 METHODS

=cut


=pod

=head1 new

$item = Genome::Model::Tools::Pcap::Item->new(children => \%children, tags => \@tags, position => $position, length => $length);
    
children - optional, a hash containing the items children, indexed by the child names.

tags - optional, an array of tags belonging to the item.

position - optional, the position of the item in the parent string.

length - optional, the length of the item in the parent item in padded base units.

=cut

sub new {
    croak("$pkg:new:no class given, quitting") if @_ < 1;
    my ($caller, %params) = @_;
    my $caller_is_obj = ref($caller);
    my $class = $caller_is_obj || $caller;
    my $self = {%params};
    bless ($self, $class);        
    if(exists $params{callbacks})
	{
		$self->callbacks($params{callbacks});
	}
    return $self;
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

sub id
{
    my ($self, $id) = @_;
    $self->{id} = $id if(@_ > 1);
    return $self->{id};
}

sub name 
{
	my ($self, $value) = @_;
	my ($name) = (caller(0))[3] =~ /.+::(.+)/;
	if(@_>1)
    {	
	   return $self->SUPER::name($value);
    }
    return $self->SUPER::name;
    #if(@_ > 1)
    #{
    #    $self->{name} = $value;
    #}
    #return $self->{name};
}

sub position 
{
	my ($self, $value) = @_;
	my ($name) = (caller(0))[3] =~ /.+::(.+)/;
	if(@_>1)
    {	
	   return $self->check_and_load_data($name, $value);
    }
    return $self->check_and_load_data($name);
}


sub length 
{
	my ($self, $value) = @_;
	my ($name) = (caller(0))[3] =~ /.+::(.+)/;
		
    if(@_>1)
    {   
       return $self->check_and_load_data($name, $value);
    }
    return $self->check_and_load_data($name);
}


sub tags
{
	my ($self, $value) = @_;
	my ($name) = (caller(0))[3] =~ /.+::(.+)/;
	
    if(@_>1)
    {   
        return $self->check_and_load_data($name, $value)||[];
    }
    return $self->check_and_load_data($name)||[];	
}

=pod

=head1 copy

$item->copy($item)           

This returns a deep copy of the $item.

=cut

sub copy
{
    my ($self,$item) = @_;
    
	return Storable::dclone($item);    
}

=pod

=head1 add_tag
    
$item->add_tag($tag)           

This adds a tag to the item's list of tags.

=cut

sub add_tag
{
    my ($self, $tag) = @_;
    
    push @{$self->tags}, $tag;
}    

=pod

=head1 copy_tag

$item->copy_tag($tag)           

This returns a copy of the tag.

=cut

sub copy_tag 
{
	my ($self,$tag) = @_;
    return Storable::dclone($tag);
}

=pod

=head1 start_position

$item->start_position           

This is a getter for the item's start position in the parent in padded units.

=cut 

sub start_position
{
	my ($self, $value) = @_;
	my ($name) = (caller(0))[3] =~ /.+::(.+)/;
    if(@_>1)
    {
        warn "This method will be read only in the future, please use the position method to set the start position\n";
       return $self->check_and_load_data('position', $value);
    }
    return $self->check_and_load_data('position');
}
=pod

=head1 add_tag

$item->end_position           

This is a getter for the item's end position in the parent in padded units.

=cut 
sub end_position 
{
    my ($self) = @_;
    
    return $self->position+$self->length;    
}

sub freeze
{
	my ($self) = @_;
	$self->SUPER::freeze;
	if($self->already_loaded("children"))
	{
		my $children = $self->children;
		foreach my $child (values %{$children})
		{
			$child->freeze;
		}
	}
}

sub thaw
{
	my ($self, $obj) = @_;
	$self->SUPER::thaw($obj);
	if($self->already_loaded("children"))
	{
		my $children = $self->children;
		foreach my $child (values %{$children})
		{
			$child->thaw($obj);
		}
	}
}

sub get_child_position_from_parent_position
{
	my ($self, $parent_pos) = @_;
	return ( $parent_pos + ( 1 - $self->position ) );    #removed -1 at the end
}

sub get_parent_position_from_child_position
{
	my ($self, $child_pos) = @_;
	return ( ( $self->position - 1 ) + $child_pos );    #removed +1 at the end

}

1;
