
# review jlolofie
#
# it is very inefficent to use this in large data sets, but
# we think the slowness is in the underlying UR code that is
# making/destroying objects while while the window is scrolling 
# around.


package Genome::Utility::Window;

use strict;
use warnings;
use Genome; 

class Genome::Utility::Window{
    is => 'UR::Object',
    has => [
        iterator => {is => 'UR::Object::Iterator'},
        range   => { is => 'Integer', default_value => 0 },

        _start => { 
            is => 'Integer', 
            default_value => 0,
            is_optional => 1,
        },
        _stop => { 
            is => 'Integer', 
            default_value => 0,
            is_optional => 1,
        },
        _max => { 
            is => 'Integer',
            is_optional => 1,
        },
        _min => { 
            is => 'Integer',
            is_optional => 1,
        },

        _iterator_position => { 
            is => 'Integer', 
            default_value => 0,
            is_optional => 1,
        },
        _iterator_done => { 
            is => 'Boolean',
            is_optional => 1, 
        },

        #_object_is_in_range => { is => 'boolean'},
        _leftover_object => { 
            is => 'UR::Object',
            is_optional => 1,
        },

        _objects => { 
            is => 'UR::Object', 
            is_many => 1,
            is_optional => 1,
        },
    ]
};

sub create
{
    my $class = shift;
    my $self = $class->SUPER::create(@_);

    $self->_start(0);
    $self->_stop(0);
    $self->_iterator_position(0);
    return $self;

=cut  
#TODO maybe replace _object_is_in_range sub w/ a closure like the following
    my $start_method = $self->object_start_method;
    my $stop_method = $self->object_stop_method;

    $self->_object_is_in_range
    (
        sub
        {
            my $object = shift;
            return -1 if $object->$stop_method < $self->_min;
            return 1 if $object->$start_method > $self->_max;
            return 0;
        }
    );
=cut

}

sub objects
{
    return @{ shift->_objects };
}

sub object_start_method
{
    return 'start';
}

sub object_stop_method
{
    return 'stop';
}

sub max
{
    return shift->_max;
}

sub min
{
    return shift->_min;
}

sub scroll
{
    my ($self, $start, $stop) = @_;

    $self->error_message("Need position to scroll") unless defined $start;

    $stop = $start unless defined $stop;
    
    return  $self->_objects  unless $self->_set_ranges($start, $stop);
    
    $self->_remove_objects;
    $self->_add_objects;
    
    return  $self->_objects ;
}

sub validate_objects
{
    my $self = shift;

    foreach my $object (  $self->_objects  )
    {
        $self->error_message("Object not in range") unless $self->_object_is_in_range($object) == 0;
        #$self->error_message("Object not in range") unless $self->_object_is_in_range->($object) == 0;
    }

    return 1;
}

sub _set_ranges
{
    my ($self, $start, $stop) = @_;

    my $current_start = $self->_start;
    my $current_stop = $self->_stop;

    return if $current_start == $start and $current_stop == $stop;

    $self->error_message("Start position to scroll to ($start) is less than the current start position ($current_start)") and die if  $start < $current_start;

    $self->error_message("Stop position to scroll to ($stop) is less than the current stop position ($current_stop)") and die if $stop < $current_stop;

    $self->_start($start);
    $self->_stop($stop);
    $self->_min( $start - $self->range );
    $self->_max( $stop + $self->range );

    return 1;
}

sub _remove_objects
{
    my $self = shift;

    my @objects;
    foreach my $object (  $self->_objects  )
    {
        push @objects, $object if $self->_object_is_in_range($object) eq 0;
        #push @objects, $object if $self->_object_is_in_range->($object) eq 0;
    }

    return $self->_objects(\@objects);
}

sub _add_objects
{
    my $self = shift;

    if ( my $leftover_object = $self->_leftover_object )
    {
        my $in_range = $self->_object_is_in_range($leftover_object); 
        #my $in_range = $self->_object_is_in_range->($leftover_object); 

        if ( $in_range == 1 ) # ahead
        {
            # done. return
            return;
        }
        elsif ( $in_range == 0 ) # in range
        {
            # add leftover to objects, contignue

            $self->_leftover_object(undef);
            $self->_objects([$self->_objects, $leftover_object]);
            #push @{ $self->_objects }, $leftover_object;
        }
        # else behind, continue
    }

    while ( 1 )
    {
        my $object = $self->_next_from_iterator
            or return;

        my $in_range = $self->_object_is_in_range($object);
        #my $in_range = $self->_object_is_in_range->($object);
        if ( $in_range == 1 ) # ahead
        {
            # done.  save object and return
            $self->_leftover_object($object);
            return;
        }
        elsif ( $in_range == 0 ) # in range
        {
            # add to objects, continue
            $self->_objects([$self->_objects, $object]);
            #push @{ $self->_objects }, $object;
        }
        # else behind, continue
    }
}

sub _next_from_iterator
{
    my $self = shift;

    return if $self->_iterator_done;
    
    my $object = $self->iterator->next;

    unless ( $object )
    {
        $self->_iterator_done(1);
        return;
    }

    $self->_iterator_position( $self->_iterator_position + 1 );

    return $object;
}

sub _object_is_in_range
{
    my ($self, $object) = @_;
    
    my $stop_method = $self->object_stop_method;
    
    return -1 if $object->$stop_method < $self->_min;
    
    my $start_method = $self->object_start_method;

    return 1 if $object->$start_method > $self->_max;

    return 0;
}

1;

=pod

=head1 Name

ModuleTemplate

=head1 Synopsis

=head1 Usage

=head1 Methods

=head2 

=over

=item I<Synopsis>

=item I<Arguments>

=item I<Returns>

=back

=head1 See Also

=head1 Disclaimer

Copyright (C) 2005 - 2008 Washington University Genome Sequencing Center

This module is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY or the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

=head1 Author(s)

B<Eddie Belter> I<ebelter@watson.wustl.edu>

=cut

#$HeadURL$
#$Id$

