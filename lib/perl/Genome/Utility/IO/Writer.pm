package Genome::Utility::IO::Writer;

#:eclark 11/17/2009 Code review.

# Shouldn't inherit from UR::Object.  Implementing new() is inconsistant with other UR based classes.
# This class should probably be the base for all file writers.

use strict;
use warnings;

use Genome;

require Cwd;
require IO::Handle;

class Genome::Utility::IO::Writer {
    is => 'UR::Object',
    is_abstract => 1,
    has_optional => [
    output => {
        type => 'String',
        doc => 'Output (file, if from command line) to write.  Defaults to STDOUT',
    },
    ],
};

sub get_original_output { # Allow getting of original output (file)
    return $_[0]->{_original_output};
}

sub new { 
    my $class = shift;
    return $class->create(@_);
}

sub create {
    my ($class, %params) = @_;

    unless ( $class->can('write_one') ) {
        $class->error_message("Can't write because there isn't a 'write_one' method in class ($class)");
        return;
    }

    my $output = delete $params{output}; # UR will complain if output doesn't match 'is' in class def

    my $self = $class->SUPER::create(%params)
        or return;
    
    unless ( $output ) { # STDOUT
        my $handle = IO::Handle->new();
        $handle->fdopen(fileno(STDOUT), "w");
        $handle->autoflush;
        $self->output($handle);
    }
    elsif ( my $output_class = ref($output) ) { # Some object?
        for my $required_method (qw/ print /) {
            unless ( $output_class->can($required_method) ) {
                $self->error_message("Output class ($output_class) can't do required method ($required_method)");
                return;
            }
        }
        $self->output($output);
    }
    else { # Assume it is a file
        $output = Cwd::abs_path($output);
        my $fh = eval { Genome::Sys->open_file_for_writing($output) };
        if (!$fh or $@) {
            $self->error_message("Can't open file $output for writing: $@");
            return;
        }
        $fh->autoflush;
        $self->{_original_output} = $output;
        $self->output($fh);
    }

    return $self;
}

sub print {
    return $_[0]->output->print($_[1]);
}

1;

=pod

=head1 Name

Genome::Utility::IO::Writer;

=head1 Synopsis

Abstract stream based writer.

=head1 Usage

B<In your class>

 package Album::Writer;

 use strict;
 use warnings;

 use above 'UR'; # or Genome

 # Declare your class UR style
 class Album::Writer {
    is => 'Genome::Utility::IO::Writer',
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

 # REQUIRED, put you writin' code here
 sub write_one {
    my ($self, $album) = @_;

    # put your print statements here, ex:
    $self->print(
        sprintf(
            "Title: %s\nArtist: %s\nGenre: %s\n\n",
            $albumn->title,
            $albumn->artist,
            $albumn->genre,
    );
    
    return 1; # return true!
 }

B<** in the code **>

 use Album::Writer;
 use IO::File;

 my $writer = Album::Writer->create(
    output => "myalbums.txt", # file or object that can 'getline' and 'seek'
 )
    or die;

 $writer->write_one(
    {
        title => 'Kool and the Gang',
        artist => 'Kool and the Gang',
        tracks => [ { name => "Celebration", rating => 5 length => '5:55' }, etc... ],
    }
 )
    or die;

=head1 Methods to Provide in Subclasses

=head2 write_one

$writer->write_one($obj)
    or die;

=over

=item I<Synopsis>   Writes the ref to the output.

=item I<Params>     ref (scalar)

=item I<Returns>    the return value of the write_one method (boolean)

=back

=head1 Methods Provided

=head2 print

 $writer->print("$string\n");

=over

=item I<Synopsis>   Prints the string to the output

=item I<Params>     none

=item I<Returns>    result of the print (boolean)

=back

=head1 See Also

I<UR>

=head1 Disclaimer

Copyright (C) 2008 Washington University Genome Sequencing Center

This module is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY or the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

=head1 Author(s)

B<Eddie Belter> <ebelter@genome.wustl.edu>

=cut

#$HeadURL$
#$Id$
