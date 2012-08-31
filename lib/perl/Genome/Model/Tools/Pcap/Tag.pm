package Genome::Model::Tools::Pcap::Tag;

use strict;
use warnings;

use base qw(Class::Accessor::Fast);

use Carp;
use Storable;

Genome::Model::Tools::Pcap::Tag->mk_accessors
(qw/
    type
    date
    scope
    source
    text
    parent
    start
    stop
    unpad_start
    unpad_stop
    no_trans
    data
    comment
    /);

our $VERSION = 0.01;

sub new
{
    croak "Genome::Model::Tools::Pcap::Tag->new no class given, quitting" unless @_;

    my ($caller, %arg) = @_;
    
    my $caller_is_obj = ref $caller;

    my $class = $caller_is_obj || $caller;

    return bless \%arg, $class;
}

sub copy
{
    my ($self, $item) = @_;
    
	return Storable::dclone($item);    
}

1;

