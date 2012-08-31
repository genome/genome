package Finfo::Reader;

use strict;
use warnings;
no warnings 'reserved';

use Finfo::Std;
require Scalar::Util;

my %io :name(io:r)
    :clo('in')
    :desc('Input file to read');
my %return_as_objs :name(return_as_objs:o) 
    :access(ro);
my %file :name(_file:p);

sub START
{
    my $self = shift;
    
    unless ( $self->can('_next') )
    {
        $self->fatal_msg("Can't read because there isn't a '_next' method");
        return;
    }
    
    if ( $self->return_as_objs and not UNIVERSAL::can($self, '_return_class') )
    {
        # TODO check return as objects data structure??
        $self->fatal_msg("Requested to return as objects, but no return class defined");
        return;
    }

    my $io = $self->io;
    if ( ref($io) )
    {
        Finfo::Validate->validate
        (
            attr => 'io for reading',
            value => $io,
            isa => 'object IO::Handle IO::String IO::File IO::Socket',
            obj => $self,
            msg => 'fatal',
        );
    }
    else
    {
        Finfo::Validate->validate
        (
            attr => 'file for reading',
            value => $io,
            isa => 'file_r',
            obj => $self,
            msg => 'fatal',
        );

        my $fh = IO::File->new("< $io");
        $self->fatal_msg("Can\'t open file ($io): $!")
            and return unless $fh;
        $self->_file($io);
        $self->io($fh);
    }
    
    return 1;
}

sub _getline : RESTRICTED
{
    return shift->io->getline;
}

sub _reset : RESTRICTED
{
    return shift->io->seek(0, 0);
}

sub next
{
    my $self = shift;

    my $ref = $self->_next;
    
    return unless defined $ref and %$ref;

    return $ref unless $self->return_as_objs;

    my $return_class = $self->_return_class;

    unless ( ref($return_class) eq 'HASH' )
    {
        return $return_class->new(%$ref);
    }

    my $class = delete $return_class->{CLASS};
    $self->error_msg("No main class given to return as objects")
        and return unless $class;
    
    foreach my $attr ( keys %$return_class )
    {
        next unless exists $ref->{attr};

        my $attr_class = $return_class->{$attr};

        my $attr_ref = ref($ref->{$attr});

        my @values = ( $attr_ref eq 'ARRAY' )
        ? @{ $ref->{$attr} }
        : $ref->{$attr};

        my @attr_objs;
        foreach my $value ( @values )
        {
            return unless Finfo::Validate->validate
            (
                attr => "Attrs for $attr_class",
                value => $ref->{$attr},
                ds => 'aryref',
                msg => "fatal_msg",
            );

            push @attr_objs, $attr_class->new(%$value)
                or return;
        }

        ( $attr_ref eq 'ARRAY' )
        ? $ref->{$attr} = \@attr_objs
        : $ref->{$attr} = $attr_objs[0];
    }

    return $class->new(%$ref);
}

sub all
{
    my $self = shift;

    my @items;
    while ( my $item = $self->next )
    {
        return if $self->short_error_msg;
        push @items, $item;
    }
    
    return if $self->short_error_msg;

    return @items;
}

1;

=pod

=head1 Name

Finfo::Reader

=head1 Synopsis

Generic base class stream based reader.  Can return values as a hashref or object.

=head1 Usage

B<In your class>

 package Album::Reader;

 use strict;
 use warnings;
no warnings 'reserved';

 use base 'Finfo::Reader';

 # use the return classes, OPTIONAL, if allowing to return_as_objects
 use Album; 
 use Album::Track; 

 # add your attributes, OPTIONAL, see Finfo::Std
 my %albums_count :name(album_count) :isa(int) default(0);
 
 # add a START method, OPTIONAL, see Fifno::Std;
 sub START # optional
 {
    my $self = shift;

    unless ( $self->_check_something )
    {
        $self->error_msg("Something went wrong");
        return;  # return false (undef or 0) to indicate error
    }
    
    return 1; # return true!!
 }

 #add a _return_class method to create a class for the parsed params, OPTIONAL
 sub _return_class 
 {
    return 'Album'; # a string
     OR
    return # a hashref
    {
        CLASS => 'Album', # main class, REQUIRED
        tracks => 'Album::Tracks # attribute (key) w/ class (value)
        etc...
    };
 }

 sub _next # required, put your parsing logic here
 {
    my $self = shift;

    my ($title, $artist, $num_of_tracks) = split(/\s+/, $self->io->getline);
    for (1..$num_of_tracks)
    {
        my $track_tokens = split(/\s+/, $self->io->getline);
        my %track_params;
        @track_params{qw/ name length rating /} = @track_tokens;

        push @tracks, \%track_params;
    }
    
    return # must be hashref
    {
        title => $title,
        artist => $artist,
        tracks => \@tracks,
    };
 }

B<in the code>

 use Album::Reader;
 use IO::File;

 my $fh = IO::File->new("< myalbmus.txt")
 my $writer = AlbumReader->new
 (
    io => $fh, 
    return_as_objects => 1
 )
    or die;

 my $album = $reader->next;

 print sprintf('%s by the famous %s', $album->title, $album->artist),"\n";

 my $track_count = 1;
 foreach my $track ( $album->tracks )
 {
    print sprintf('#%d %s rated %s', $track_count++; $track->name, $track->rating),"\n";
 }
 # etc...

=head1 Methods to Overwrite in Your Class

=head2 _next

=over

=item I<Synopsis>   This is where the particular parsing code goes.

=item I<Params>     none

=item I<Returns>    hashref

=back

=head2 _return_class

=over

=item I<Synopsis>   This is the return class of the objects.  If defined and return_as_objs is true, the next method will automatically return the hashref as the object.

=item I<Params>     none

=item I<Returns>    string OR hashref (see above)

=back

=head1 Public Base Class Methods - Do not overwrite!

=head2 next

my $ref(or object) = $reader->next;

=over

=item I<Synopsis>   Gets the next ref/object form the io.  Calls _next in your class.

=item I<Params>     none

=item I<Returns>    scalar (hashref or object)

=back

=head2 all

my @refs (or objects) = $reader->all;

=over

=item I<Synopsis>   Gets all the refs/objects form the io.  Calls _next in your class until it returns undefined or an error is encountered

=item I<Params>     none

=item I<Returns>    array (hashrefs or objects)

=back

=head1 Private Methods Restricted to Child Classes

=head2 _getline

$reader->_getline
    or die;

=over

=item I<Synopsis>   Gets the next line form the io

=item I<Params>     none

=item I<Returns>    the line

=back

=head2 _reset

$reader->_reset
    or die;

=over

=item I<Synopsis>   Resets (seek) the io to the beginning

=item I<Params>     none

=item I<Returns>    the result of the io->seek

=back

=head1 See Also

I<Finfo::Std>, I<Finfo::Logging>

=head1 Disclaimer

Copyright (C) 2007 Washington University Genome Sequencing Center

This module is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY or the implied warranty of MERCHANTABILITY
or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public
License for more details.

=head1 Author(s)

Eddie Belter <ebelter@watson.wustl.edu>

=cut

#$HeadURL$
#$Id$
