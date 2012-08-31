package Genome::Model::Tools::Pcap::Sources::Contig;
our $VERSION = 0.01;

use strict;
use warnings;
use Carp;
use Storable;
use base (qw(Genome::Model::Tools::Pcap::Sources::SequenceItem));


#Contig Data Structure
#Contig:
#	sequence (bases and quality)
#	children 
#	get_child
#	padded_base_count
#	base_count
#	read_count
#	complemented
#	tags


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
	my $name = (caller(0))[3];
    croak "$name is an abstract base method!\n";
}

sub base_count
{
	my $name = (caller(0))[3];
    croak "$name is an abstract base method!\n";
}

sub base_segment_count
{
	my $name = (caller(0))[3];
    croak "$name is an abstract base method!\n";
}

sub complemented
{
	my $name = (caller(0))[3];
    croak "$name is an abstract base method!\n";
}

1;
