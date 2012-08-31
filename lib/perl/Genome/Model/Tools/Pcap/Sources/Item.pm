package Genome::Model::Tools::Pcap::Sources::Item;
our $VERSION = 0.01;

use strict;
use warnings;
use Carp;

sub new 
{
	my $name = (caller(0))[3];
    croak("__PACKAGE__:new:no class given, quitting") if @_ < 1;
	my ($caller, %params) = @_; 
    my $caller_is_obj = ref($caller);
    my $class = $caller_is_obj || $caller;      
	my $self = \%params; 
	bless ($self, $class);
	return $self;
}

#TODO: ace specific, the two methods below need to be moved to Genome::Model::Tools::Pcap::Sources::Ace
sub freeze
{
	
}

sub thaw
{
	
}

sub get_map {
	my $name = (caller(0))[3];
    croak "$name is an abstract base method!\n";	
}

sub children 
{
	my $name = (caller(0))[3];
    croak "$name is an abstract base method!\n";	
}

sub name 
{
	my $name = (caller(0))[3];
    croak "$name is an abstract base method!\n";
}

sub position 
{
	my $name = (caller(0))[3];
    croak "$name is an abstract base method!\n";	
}


sub length 
{
	my $name = (caller(0))[3];
    croak "$name is an abstract base method!\n";
}

sub tags
{
	my $name = (caller(0))[3];
    croak "$name is an abstract base method!\n";
}

sub copy
{
	my $name = (caller(0))[3];
    croak "$name is an abstract base method!\n";
}

sub add_tag
{
	my $name = (caller(0))[3];
    croak "$name is an abstract base method!\n";
} 

sub copy_tag 
{
	my $name = (caller(0))[3];
    croak "$name is an abstract base method!\n";
}

sub start_position
{
	my $name = (caller(0))[3];
    croak "$name is an abstract base method!\n";
}

sub end_position 
{
	my $name = (caller(0))[3];
    croak "$name is an abstract base method!\n";
}

1;
