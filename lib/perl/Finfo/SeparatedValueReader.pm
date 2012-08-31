package Finfo::SeparatedValueReader;

use strict;
use warnings;
no warnings 'reserved';

use base 'Finfo::Reader';

use Data::Dumper;

my %headers :name(headers:o)
    :ds(aryref) 
    :isa(string)
    :desc('Headers for the file.  If none given, they will be retrieved from the file.');
my %sep :name(separator:o)
    :isa(string)
    :default(',')
    :desc('The value separator character');
my %is_regex :name(is_regex:o)
    :isa(boolean)
    :default(0)
    :desc('interprets separator as regex'); 
my %ln :name(_line_number:p)
    :isa('int')
    :default(0);
my %split :name(_split:p)
    :isa('code');


sub START {
    my $self = shift;
    
    my $sep = $self->separator;
    if ($self->is_regex){ #Adding -1 as the LIMIT argument to split ensures that the correct # of values on the line are returned, regardless of empty trailing results
        $self->_split(sub{my $string = shift; return split(/$sep/, $string, -1)});
    }else{
        $self->_split(sub{my $string = shift; return split(/\Q$sep\E/, $string, -1)});
    }
    return 1 if $self->headers;

    my @headers = $self->_get_and_split_next_line;
    $self->error_msg("No headers found in io")
        and return unless @headers;
    
    return $self->headers(\@headers);
}

sub line_number {
    return shift->_line_number;
}

sub _next {
    my $self = shift;

    my @values = $self->_get_and_split_next_line
        or return;
    foreach (@values){
        $_ =~ s/^\s*['"]?//;
        $_ =~ s/['"]?\s*$//;
    }

    my $headers = $self->headers;
    
    $self->fatal_msg (
        sprintf(
            'Expected %d values, got %d on line %d in %s', 
            scalar @$headers, 
            scalar @values,
            $self->_line_number,
            ( $self->_file || ref $self->io ),
        )
    ) unless scalar @$headers == scalar @values;
    
    my %data;
    @data{@$headers} = @values;
    #TODO remove quotes? 
    #@data{@$headers} =
    
    return \%data;
}

sub _get_and_split_next_line {
    my $self = shift;

    my $line = $self->_getline
        or return;

    $self->_line_number( $self->_line_number + 1 );
    
    chomp $line;

    return $self->_split->($line); 
}

1;

=pod

=head1 Name

Finfo::SeparatedValueReader

=head1 Synopsis

Iterates through an io, splitting each line by the given separator.  If no headers are given, they will be derived from the first line of the io, being split by the separator.

=head1 Usage

 use Finfo::SeparatedValueReader;

 my $reader = Finfo::SeparatedValueReader->new (
    io => 'myalbmus.txt', # req; file or IO:: object
    headers => [qw/ title artist /], # opt; headers for the file
    separator => '\t', # opt; default is ','
    is_regex => 'set this flag if your separator is a regular expression, otherwise the literal characters of the separator will be used'
 );

 while ( my $album = $reader->next ) {
    print sprintf('%s by the famous %s', $album->{title}, $album->{artist}),"\n";
 }

=head1 Methods

=head2 next

 my $ref = $reader->next;

=over

=item I<Synopsis>   Gets the next ref form the io.

=item I<Params>     none

=item I<Returns>    scalar (hashref)

=back

=head2 all

 my @refs = $reader->all;

=over

=item I<Synopsis>   Gets all the objects form the io.

=item I<Params>     none

=item I<Returns>    array (hashrefs)

=back

=head1 See Also

I<Finfo::Reader> (inherits from), I<Finfo::Std>, I<Finfo::Logging>

=head1 Disclaimer

Copyright (C) 2008 Washington University Genome Sequencing Center

This module is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY or the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

=head1 Author(s)

Eddie Belter <ebelter@watson.wustl.edu>

=cut

#$HeadURL$
#$Id$
