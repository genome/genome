package Finishing::Assembly::Iterator;

use strict;
use warnings;
no warnings 'reserved';

use Finfo::Std;

use Data::Dumper;

my %constructor 
    :name(object_constructor:r);
my %iterator
    :name(iterator:r)
    :isa('object DBIx::Class::ResultSet Finfo::Iterator');

sub next
{
	my $self = shift;

    my $obj = $self->iterator->next;

    return unless $obj;

    return $self->object_constructor->($obj);
}

sub count 
{
    my $self = shift;

    return $self->iterator->count(@_);
}

sub all 
{
	my $self = shift;
    
    my @objects = $self->iterator->all;

    return unless @objects;
    
    return map { $self->object_constructor->($_) } @objects;
}

sub reset 
{
	my $self = shift;

    return $self->iterator->reset;
}

sub first 
{
	my $self = shift;

    my $obj = $self->iterator->first;
    
    return unless $obj;

    return $self->object_constructor->($obj);
}

sub find
{
	my $self = shift;

    my $obj = $self->iterator->find(@_);
    
    return unless $obj;

    return $self->object_constructor->($obj);
}

sub search
{
    my $self = shift;

    my $iterator = $self->iterator->search(@_);

    if ( wantarray ) 
    {
        my @objects = $iterator->all;
        return unless @objects;
        return map { $self->object_constructor->($_) } @objects;
    }
    else
    {
        return __PACKAGE__->new
        (
            object_constructor => $self->object_constructor,
            iterator => $iterator,
        );
    }
}

1;

#$HeadURL$
#$Id$
