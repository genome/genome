package Genome::Utility::IO::SeparatedValueReader;

use strict;
use warnings;

use Genome;

class Genome::Utility::IO::SeparatedValueReader {
    is => 'Genome::Utility::IO::Reader', 
    has_optional => [
        headers => {
            type => 'Array',
            doc => 'Headers for the file.  If none given, they will be retrieved from the input.'
        },
        separator => {
            type => 'String',
            default => ',',
            doc => 'The value of the separator character.  Default: ","'
        },
        is_regex => {
            type => 'Boolean',
            default => 0,
            doc => 'Interprets separator as regex'
        },
        ignore_extra_columns => {
            type => 'Boolean',
            default => 0,
            doc => 'Rather than crash when extra columns are found, just drop/ignore them.'
        },
        ignore_lines_starting_with => {
            type => 'Text',
            doc => 'Ignore initial lines starting with character(s).'
        },
        current_line => {
            type => 'Text',
            doc => 'The current original line of input from the file, pre splitting.'
        },
    ],
};

sub line_number {
    return $_[0]->{_line_number};
}

sub reset { 
    my $self = shift;

    $self->SUPER::reset;
    $self->{_line_number} = 0;

    if ( $self->_headers_were_in_input ) { # skip to data
        $self->getline;
    }

    return 1;
}

sub getline { 
    my $self = shift;

    my $line = $self->SUPER::getline
        or return;

    $self->_increment_line_number;

    return $line;
}

sub create {
    my $class = shift;
    my %params = @_;

    my $headers = delete $params{headers}; # prevent UR from sorting our headers!
    my $self = $class->SUPER::create(%params)
        or return;

    my $sep = $self->separator;

    my $regexp;
    if ($self->is_regex){ 
        # Adding -1 as the LIMIT argument to split ensures that the correct # of values on the line 
        #  are returned, regardless of empty trailing results
        $regexp = qr/$sep/;
    }
    else {
        $regexp = qr/\Q$sep\E/;
    }
    $self->{_split} = sub{ return _strip_quotes(split($regexp, $_[0], -1)) };

    if (defined($self->ignore_lines_starting_with)) {
        my $char = $self->ignore_lines_starting_with;
        my $line = $self->getline;
        my $offset = length($line);
        while ($line =~ /^$char/) {
            $line = $self->getline;
            $offset = length($line);
        }
        my $begin = tell($self->input) - $offset;
        $self->input->seek($begin,0);
    }
    
    if ( $headers ) {
        $self->headers($headers);
    }
    else {
        my $line = $self->getline;
        my @headers = $self->_splitline($line);
        $self->error_message("No headers found in io")
            and return unless @headers;
        $self->headers(\@headers);
        $self->{_headers_were_in_input} = 1;
    }

    return $self;
}

sub next {
    my $self = shift;
    my $line = $self->getline or return;
    my @values = $self->_splitline($line)
        or return;

    unless ( @{$self->headers} == @values ) {

        # If we dont care about extra columns, all is well... unless we dont at least have the minimum required
        if (!($self->ignore_extra_columns)||(@{$self->headers} > @values)) {
            #  Bomb out if we dont want extra columns
            $self->error_message(
               sprintf(
                    'Expected %d values, got %d on line %d in %s', 
                    scalar @{$self->headers}, 
                    scalar @values,
                    $self->line_number,
                    ( $self->get_original_input || ref $self->io ),
                )
            );
            return;
        }
    }
    
    my %data;
    @data{ @{$self->headers} } = @values;
    $self->current_line($line);
    return \%data;
}

sub _splitline {
    my $self = shift;
    my $line = shift;
    chomp $line;

    return $self->{_split}->($line); 
}

sub _increment_line_number {
    return $_[0]->{_line_number}++;
}

sub _headers_were_in_input {
    return $_[0]->{_headers_were_in_input};
}

my $strip_leading_quotes_regex = qr/^\s*['"]?/;
my $strip_trailing_quotes_regex = qr/['"]?\s*$/;
sub _strip_quotes {
    map { $_ =~ s/$strip_leading_quotes_regex//; $_ =~ s/$strip_trailing_quotes_regex//; $_ } @_;
}


1;

=pod

=head1 Name

Genome::Utility::IO::SeparatedValueReader

=head1 Synopsis

A stream based reader that splits each line by the given separator.  If no headers are given, they will be derived from the first line of the io, being split by the separator.

=head1 Usage

 use Genome::Utility::IO::SeparatedValueReader;

 my $reader = Genome::Utility::IO::SeparatedValueReader->new (
    input => 'albums.txt', # REQ: file or object that can 'getline' and 'seek'
    headers => [qw/ title artist /], # OPT; headers for the file
    separator => '\t', # OPT; default is ','
    is_regex => 1, # OPT: 'set this flag if your separator is a regular expression, otherwise the literal characters of the separator will be used'
 );

 while ( my $album = $reader->next ) {
    print sprintf('%s by the famous %s', $album->{title}, $album->{artist}),"\n";
 }

=head1 Methods 

=head2 next

 my $ref = $reader->next;

=over

=item I<Synopsis>   Gets the next hashref form the input.

=item I<Params>     none

=item I<Returns>    scalar (hashref)

=back

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

=head2 line_number

 my $line_number = $reader->line_number;

=over

=item I<Synopsis>   Gets the current line number (position) of the input

=item I<Params>     none

=item I<Returns>    line numeber (int)

=back

=head1 See Also

I<Genome::Utility::IO::SeparatedValueReader> (inherits from), I<UR>, I<Genome>

=head1 Disclaimer

Copyright (C) 2008 Washington University Genome Sequencing Center

This module is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY or the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

=head1 Author(s)

Eddie Belter <ebelter@watson.wustl.edu>

=cut

