package Genome::Model::Tools::Pcap::Sources::SequenceItem;
our $VERSION = 0.01;

use strict;

use warnings;
use Carp;
use Storable;
use Genome::Model::Tools::Pcap::Transform;
use base(qw(Genome::Model::Tools::Pcap::Sources::Item));

sub new {
    croak("__PACKAGE__:new:no class given, quitting") if @_ < 1;
    my ($caller, %args) = @_;
    my $caller_is_obj = ref($caller);
    my $class = $caller_is_obj || $caller;
    my $self = \%args;
    bless ($self, $class);		
	     
    return $self;
}

sub get_map {
	my $name = (caller(0))[3];
    croak "$name is an abstract base method!\n";
}

sub _transform
{
	my $name = (caller(0))[3];
    croak "$name is an abstract base method!\n";
}

sub get_transform
{
	my $name = (caller(0))[3];
    croak "$name is an abstract base method!\n";
}

sub _load_transform
{
	my $name = (caller(0))[3];
    croak "$name is an abstract base method!\n";
}

sub _has_alignment
{
	my $name = (caller(0))[3];
    croak "$name is an abstract base method!\n";
}

sub padded_base_string
{
	my $name = (caller(0))[3];
    croak "$name is an abstract base method!\n";
}

sub padded_base_quality
{
	my $name = (caller(0))[3];
    croak "$name is an abstract base method!\n";
}

sub unpadded_base_string
{
	my $name = (caller(0))[3];
    croak "$name is an abstract base method!\n";
}

sub unpadded_base_quality
{
	my $name = (caller(0))[3];
    croak "$name is an abstract base method!\n";
}

sub get_padded_base_quality
{
	my $name = (caller(0))[3];
    croak "$name is an abstract base method!\n";
}

sub get_padded_base_value
{
	my $name = (caller(0))[3];
    croak "$name is an abstract base method!\n";
}

sub has_alignment
{
	my $name = (caller(0))[3];
    croak "$name is an abstract base method!\n";
}

sub copy
{
    my ($self,$item) = @_;
    
	return Storable::dclone($item);    
}

sub length
{
	my $name = (caller(0))[3];
    croak "$name is an abstract base method!\n";
}

1;
