package Genome::Model::Tools::Pcap::Sources::Ace::SequenceItem;
our $VERSION = 0.01;

use strict;
use warnings;
use Carp;

use Storable;
use base (qw(Genome::Model::Tools::Pcap::Sources::SequenceItem));

sub new 
{
    croak("__PACKAGE__:new:no class given, quitting") if @_ < 1;
    my ($caller, %params) = @_; 
    my $caller_is_obj = ref($caller);
    my $class = $caller_is_obj || $caller;
    my $self = $class->SUPER::new(%params);

    $self->_index($params{index});
    return $self;
}

#methods inherited from Item Source

sub freeze
{
	my ($self) = @_;
	$self->{fh} = undef;
	$self->{reader}->{'input'} = undef;
}

sub thaw
{
	my ($self, $obj) = @_;
    return unless (defined $self->_index);
    $self->{fh} = Genome::Model::Tools::Pcap::FHManager->get_fh($self->_index()->{file_name});
	$self->{reader}->{'input'} = $self->{fh};
}

1;
