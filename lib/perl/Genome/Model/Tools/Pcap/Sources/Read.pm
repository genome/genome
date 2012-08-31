package Genome::Model::Tools::Pcap::Sources::Read;
our $VERSION = 0.01;

use strict;
use warnings;
use Carp;

use Storable;
use base (qw(Genome::Model::Tools::Pcap::Sources::SequenceItem));

=head1 NAME

Read - Read  Object.

=cut
#Read Data Structure
#Read:
#	sequence (bases and quality)
#	align_clip_start
#	align_clip_end
#	qual_clip_start
#	qual_clip_end
#	padded_base_count
#	base_count
#	complemented
#	tags
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
	my $name = (caller(0))[3];
    croak "$name is an abstract base method!\n";
}

sub complemented
{
	my $name = (caller(0))[3];
    croak "$name is an abstract base method!\n";
}

sub align_clip_start 
{
	my $name = (caller(0))[3];
    croak "$name is an abstract base method!\n";
}

sub align_clip_end 
{
	my $name = (caller(0))[3];
    croak "$name is an abstract base method!\n";
}

sub qual_clip_start 
{
	my $name = (caller(0))[3];
    croak "$name is an abstract base method!\n";
}

sub qual_clip_end 
{
	my $name = (caller(0))[3];
    croak "$name is an abstract base method!\n";
}

sub info_count 
{
	my $name = (caller(0))[3];
    croak "$name is an abstract base method!\n";
}

sub chromat_file 
{
	my $name = (caller(0))[3];
    croak "$name is an abstract base method!\n";
}

sub phd_file 
{
	my $name = (caller(0))[3];
    croak "$name is an abstract base method!\n";
}

sub time
{
	my $name = (caller(0))[3];
    croak "$name is an abstract base method!\n";
}

1;
