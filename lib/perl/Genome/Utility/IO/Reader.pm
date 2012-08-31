package Genome::Utility::IO::Reader;

#:eclark 11/17/2009 Code review.

# Shouldn't inherit from UR::Object.  Implementing new() is inconsistant with other UR based classes.
# This class should probably be the base for all file type readers.  Should be is_abstract

use strict;
use warnings;

use Genome;

class Genome::Utility::IO::Reader {
    is => 'UR::Object',
    has => [
    input => {
        type => 'String',
        is_optional => 0,
        doc => 'Input (file, if from command line) to read',
    },
    ],
};

sub get_original_input { # Allow getting of original input (file)
    return shift->{_original_input};
}

sub new { 
    my $class = shift;
    return $class->create(@_);
}

sub create {
    my ($class, %params) = @_;

    my $input = delete $params{input};
    my $self = $class->SUPER::create(%params)
        or return;

    $self->error_message("Input is required for class ($class)")
        and return unless defined $input;

    if ( my $input_class = ref($input) ) {
        for my $required_method (qw/ getline seek /) {
            unless ( $input_class->can($required_method) ) {
                $self->error_message("Input class ($input_class) can't do required method ($required_method)");
                return;
            }
        }
        $self->input($input);
    }
    else {
        my $fh = eval { Genome::Sys->open_file_for_reading($input) };
        if (!$fh or $@) {
            $self->error_message("Can't open file $input for reading: $@");
            return;
        }

        $self->{_original_input} = $input;
        $self->input($fh);
    }

    return $self;
}

sub getline {
    return shift->input->getline;
}

sub reset {
    return shift->input->seek(0, 0);
}

sub all {
    my $self = shift;

    my @items;
    while ( 1 ) {
        my $item = $self->next;
        unless  ( $item ) {
            return if $self->error_message; # had a problem
            last; # got 'em!
        }
        push @items, $item;
    }

    return @items;
}

1;

=pod

=head1 Name

Genome::Utility::IO::Reader

=head1 Synopsis

Abstract stream based reader.

=head1 Usage

B<In your class>

 package Album::Reader;

 use strict;
 use warnings;

 use above 'UR'; # or Genome

 # Declare your class UR style
 class Album::Reader {
    is => 'Genome::Utility::IO::Reader',
 };
 
 # OPTIONAL - add a create method
 sub create {
    my $class = shift;

    my $self = $class->SUPER::create(@_); # REQUIRED to call SUPER::create

    unless ( $self->_check_something ) {
        $self->error_message("Something went wrong");
        return; # return false (undef or 0) to indicate error
    }
    
    return $self;
 }

 sub next { # REQUIRED, put your parsing logic here
    my $self = shift;

    my ($title, $artist, $num_of_tracks) = split(/\s+/, $self->getline);
    for (1..$num_of_tracks) {
        my $track_tokens = split(/\s+/, $self->getline);
        my %track_params;
        @track_params{qw/ name length rating /} = @track_tokens;

        push @tracks, \%track_params;
    }
    
    return { # must be hashref
        title => $title,
        artist => $artist,
        tracks => \@tracks,
    };
 }

B<in the code>

 use Album::Reader;

 my $reader = Album::Reader->create(
    input => 'albums.txt', # file or object that can 'getline' and 'seek'
 )
    or die;

 my $album = $reader->next;

 print sprintf('%s by the famous %s', $album->title, $album->artist),"\n";

 my $track_count = 1;
 foreach my $track ( $album->tracks ) {
    print sprintf('#%d %s rated %s', $track_count++; $track->name, $track->rating),"\n";
 }
 # etc...

=head1 Methods to Provide in Subclasses

=head2 next

 my $ref (or object) = $reader->next;

=over

=item I<Synopsis>   Gets the next ref/object form the input.

=item I<Params>     none

=item I<Returns>    scalar (hashref or object)

=back

=head1 Methods Provided

=head2 all

 my @refs (or objects) = $reader->all;

=over

=item I<Synopsis>   Gets all the refs/objects form the input.  Calls _next in your class until it returns undefined or an error is encountered

=item I<Params>     none

=item I<Returns>    array (hashrefs or objects)

=back

=head2 getline

 $reader->getline
    or die;

=over

=item I<Synopsis>   Returns the next line form the input (not chomped)

=item I<Params>     none

=item I<Returns>    scalar (string)

=back

=head2 reset

 $reader->reset
    or die;

=over

=item I<Synopsis>   Resets (seek) the input to the beginning

=item I<Params>     none

=item I<Returns>    the result of the $self->input->seek (boolean)

=back

=head1 See Also

I<UR::Object>

=head1 Disclaimer

Copyright (C) 2008 Washington University Genome Sequencing Center

This module is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY or the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

=head1 Author(s)

B<Eddie Belter> <ebelter@genome.wustl.edu>

=cut

#$HeadURL$
#$Id$
