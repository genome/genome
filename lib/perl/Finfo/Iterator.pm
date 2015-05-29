package Finfo::Iterator;

use strict;
use warnings;
no warnings 'reserved';

use Finfo::Std;

use Data::Dumper;

my %objects :name(objects:o) :ds(aryref) :empty_ok(1) :access(ro);
my %ids :name(ids:o) :ds(aryref) :empty_ok(1);
my %cb :name(cb:o) :isa(code);
my %id_method :name(id_method:o) :isa(string) :default(id);
my %count :name(_count:p) :isa('int gte 0');
my %pos :name(_position:p) :isa('int gte -1') :default(-1);

sub START
{
    my $self = shift;

    my $ids = $self->ids;
    my $objects = $self->objects;

    $self->fatal_msg("Can't create iterator with both ids and objects") if $ids and $objects;
    $self->fatal_msg("Need objects or ids to create iterator") unless $ids or $objects;

    if ( $objects )
    {
        my $i = 0;
        $self->ids([ map { $i++ } @{ $self->objects } ]);
        $self->cb( sub{ $self->_get_object(@_); });
    }


    return 1;
}

sub _get_object : PRIVATE
{
    my ($self, $id) = @_;

    return $self->objects->[$id];
}

sub count
{
    my $self = shift;

    return scalar( @{ $self->ids } );
}

sub reset
{
    my $self = shift;
    
    return $self->_position(-1);
}

sub first 
{
    my $self = shift;

    return $self->cb->( $self->ids->[0] );
}

sub last 
{
    my $self = shift;
    
    return $self->cb->( $self->ids->[ ($self->count - 1) ] );
}

sub next
{
    my $self = shift;

    my $position = $self->_increment_position;
    return if $position > ($self->count - 1);

    return $self->cb->( $self->ids->[$position] );
}

sub all 
{
    my $self = shift;

    return map { $self->cb->($_) } @{ $self->ids };
}

sub find
{
    my ($self, $cond) = @_;
    
    my ($method, $val) = %$cond;
    
    if ( $method eq $self->id_method and not $self->objects )
    {
        return $self->cb->($val);
    }

    my @objects = $self->all;
    foreach my $obj ( @objects )
    {
        return $obj if $obj->$method eq $val;
    }

    return;
}

sub search
{
    my ($self, $conditions, $params) = @_;

    unless ( $conditions )
    {
        if ( wantarray )
        {
            return ( $self->objects )
            ? @{ $self->objects }
            : $self->all;
        }
        else
        {
            my $method = $self->id_method;
            if ( $self->objects )
            {
                return __PACKAGE__->new(objects => $self->objects);
            }
            else
            {
                return __PACKAGE__->new 
                (
                    ids => $self->ids, 
                    cb => $self->cb,
                    id_method => $self->id_method,
                );
            }
        }
    }

    Finfo::Validate->validate
    (
        attr => 'search conditions',
        value => $conditions,
        ds => 'hashref',
        msg => 'fatal',
    );

    if ( $params )
    {
        Finfo::Validate->validate
        (
            attr => 'search params',
            value => $params,
            ds => 'hashref',
            msg => 'fatal',
        );
    }

    my @objects = $self->all;
    my @matching_objects;
    OBJECT: foreach my $obj ( @objects )
    {
        ATTR: foreach my $attr ( keys %$conditions )
        {
            my $val = $conditions->{$attr};
            if ( exists $params->{$attr} )
            {
                my $match = $params->{$attr};
                if ( $match eq 'like' )
                {
                    next OBJECT unless $obj->$attr =~ /$val/i;
                }
                elsif ( $match eq 'gte' )
                {
                    next OBJECT unless $obj->$attr >= $val;
                }
                elsif ( $match eq 'gt' )
                {
                    next OBJECT unless $obj->$attr > $val;
                }
                elsif ( $match eq 'lte' )
                {
                    next OBJECT unless $obj->$attr <= $val;
                }
                elsif ( $match eq 'lt' )
                {
                    next OBJECT unless $obj->$attr < $val;
                }
                #TODO add more
                else
                {
                    $self->fatal_msg("Unknown match type ($match)");
                }
            }
            else
            {
                if ( ref($val) )
                {
                    my $attr_val = $obj->$attr;
                    next OBJECT unless grep { $attr_val eq $_ } @$val;
                }
                else
                {
                    next OBJECT unless $obj->$attr eq $val;
                }
            }
            push @matching_objects, $obj;
        }
    }

    if ( wantarray )
    {
        return @matching_objects;
    }
    else
    {
        my $method = $self->id_method;
        __PACKAGE__->new 
        (
            ids => [ map { $_->$method } @matching_objects ], 
            cb => $self->cb,
            id_method => $method,
        );
    }
}

sub _increment_position : PRIVATE
{
    my $self = shift;

    return $self->_position( $self->_position + 1 );
}

1;

=pod

=head1 Name

Finfo::Iterator

=head1 Synopsis

=head1 Usage


 use Finfo::Iterator;

B<create with ids and a callback to get the object associated with the id:>

 my $iterator = Finfo::Iterator->new
 (
    ids => [qw/ 1 2 3 /],
    callback => sub{ return _get_object_by_id(@_); },
 );

B<or use the objects themselves:>

 my $iterator = Finfo::Iterator->new
 (
    objects => \@objects,
 );

B<then...>

 while ( my $obj = $iterator->next )
 {
    ...
 }

=head1 Methods

=head2 count

 my $count = $iterator->count;
 
=over

=item I<Synopsis>   Gets the number of objects

=item I<Arguments>  none

=item I<Returns>    count (integer)

=back

=head2 next

 my $object = $iterator->next;

 or 

 while ( my $object = $iterator->next )
 {
    ...
 }
 
=over

=item I<Synopsis>   gets the next object, use inb a while loop

=item I<Arguments>  none

=item I<Returns>    object

=back

=head2 all

 my @objects = $iterator->all;
 
=over

=item I<Synopsis>   gets all of the objects

=item I<Arguments>  none

=item I<Returns>    objects (array)

=back

=head2 first

 my $object = $iterator->first;
 
=over

=item I<Synopsis>   gets the first object in the series

=item I<Arguments>  none

=item I<Returns>    object

=back

=head2 last

 my $object = $iterator->last;

=over

=item I<Synopsis>   gets the last object in the series

=item I<Arguments>  none

=item I<Returns>    object

=back

=head2 reset

 $iterator->reset;

=over

=item I<Synopsis>   resets the iterator to the first position

=item I<Arguments>  none

=item I<Returns>    true on success

=back

=head2 search

 my $new_iterator = $iterator->search
 (

    { name => 'tom' }, # conditions
    { name => 'like' }, # parameters
 );

 or
 
 my @objects = $iterator->search
 (
    { name => 'tom' }, # conditions
    { name => 'like' }, # parameters
 );

=over

=item I<Synopsis>   goes through the objects and determines if the meet the conditions and parameters

=item I<Arguments>  conditions (hashref of object attributes and values), params (hashref of object attributes and comparison function)

=item I<Returns>    objects if wantarry, otherwise an array of the objects

=item I<TODO>       document functions, add more functions

=back

=head2 find

 my $object = $iterator->find({ name => 'Tom' });

=over

=item I<Synopsis>   given an attribute and value, finds the first object that matches

=item I<Arguments>  conditions (hashref of an attribute and a value)

=item I<Returns>    object or undef if not found

=back

=head1 Disclaimer

Copyright (C) 2005 - 2008 Washington University Genome Sequencing Center

This module is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY or the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

=head1 Author(s)

B<Eddie Belter> I<ebelter@watson.wustl.edu>

=cut

#$HeadURL$
#$Id$
