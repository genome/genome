package Finfo::Writer;

use strict;
use warnings;
no warnings 'reserved';

use Finfo::Std;
require Scalar::Util;

my %io :name(io:r)
    :clo('out=s')
    :desc('Output file to write');
my %file :name(_file:p);

sub START
{
    my $self = shift;

    unless ( $self->can('_write_one') )
    {
        $self->fatal_msg("Can't write because there isn't a '_write_one' method");
        return;
    }
    
    my $io = $self->io;
    if ( ref($io) )
    {
        Finfo::Validate->validate
        (
            attr => 'io for writing',
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
            attr => 'file for writing',
            value => $io,
            isa => 'file_w',
            obj => $self,
            msg => 'fatal',
        );

        my $fh = IO::File->new("> $io");
        $self->fatal_msg("Can\'t open file ($io): $!")
            and return unless $fh;
        $self->_file($io);
        $self->io($fh);
    }
 
    return 1;
}

sub write_one
{
    my ($self, $ref) = @_;

    return $self->_write_one($ref);
}

sub write_many
{
    my ($self, $refs) = @_;

    Finfo::Validate->validate
    (
        attr => 'objects to write',
        value => $refs,
        ds => 'aryref',
        msg => 'fatal',
    );
    
    my $count = 0;
    foreach my $ref ( @$refs )
    {
        return unless $self->_write_one($ref);
        $count++;
    }

    return $count;
}

1;

=pod

=head1 Name

Finfo::Writer

=head1 Synopsis

Generic base class stream based writer.

=head1 Usage

B<In your class>

 package Album::Writer;

 use strict;
 use warnings;
no warnings 'reserved';

 use base 'Finfo::Writer';

 # additional attributes, OPTIONAL, see Finfo::Std
 my %album_count :name(album_count) :isa(int) :default(0);

 # add a START (or BUILD) method, OPTIONAL, see Finfo::Std;
 sub START 
 {
    my $self = shift;

    unless ( $self->_check_something )
    {
        $self->error_msg("Something went wrong");
        return;  # return false (undef or 0) to indicate error
    }
    
    return 1; # return true!!
 }

 sub _write_one # required, put you writin' code here
 {
    my ($self, $album) = @_;

    Finfo::Validate->validate # optional
    (
        attr => 'album to write',
        value => $album,
        isa => 'object Album',
    );
    
    # put your print statements here, ex:
    $self->io->print
    (
        sprintf
        (
            "Title: %s\nArtist: %s\nGenre: %s\n\n",
            $albumn->title,
            $albumn->artist,
            $albumn->genre,
    );
    
    return 1; # return true!
 }

**** in the code ****

 use Album;
 use Album::Writer;
 use IO::File;

 my $writer = Album::Writer->new
 (
    io => IO::String->new("< myalbums.txt")
 )
    or die;

 my $album = Album->new
 (
    title => 'Kool and the Gang',
    artist => 'Kool and the Gang',
    tracks => [ { name => "Celebration", rating => 5 length => '5:55' }, etc... ],
 );

 $writer->write_one($album)
    or die;

=head1 Methods - Base Class!  Overwrite in your class!

=head2 write_one

$writer->write_one($obj)
    or die;

=over

=item I<Synopsis>   Writes the ref/object to the io.  Calls _write_one in your class.

=item I<Params>     ref/object (scalar)

=item I<Returns>    the return value of the _write_one method (true on success)

=back

=head2 write_many

$writer->write_many(@objs)
    or die;

=over

=item I<Synopsis>   Writes the refs/objects.  Calls _write_one in your class foreach ref/object

=item I<Params>     hashrefs or objects (array)

=item I<Returns>    num of objects written (int)

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
