package Genome::Utility::IO::SeparatedValueWriter;

use strict;
use warnings;

use Genome;

use Data::Compare 'Compare';
use Data::Dumper 'Dumper';
require Scalar::Util;

class Genome::Utility::IO::SeparatedValueWriter {
    is => 'Genome::Utility::IO::Writer', 
    has => [
        headers => {
            type => 'Array',
            doc => 'Headers to write and use in ordering value order.'
        },
    ],
    has_optional => [
        separator => {
            type => 'String',
            default => ',',
            doc => 'The value of the separator character.  Default: ","'
        },
        in_place_of_null_value => {
            doc => 'Use this in place of an undefined value.',
        },
        print_headers => {
            is => 'Boolean',
            doc => 'A flag to print the headers on the first line.',
            default_value => 1,
        },
        ignore_extra_columns => {
            is => 'Boolean',
            doc => 'A flag to ignore extra data key/value pairs not in headers',
            default_value => 0,
        },
    ],
};

sub create {
    my ($class, %params) = @_;

    my $headers = delete $params{headers}; # prevent UR from sorting our headers!
    unless ( $headers ) {
        $class->error_message("Headers are required to create.");
        return;
    }

    my $self = $class->SUPER::create(%params)
        or return;

    $self->headers($headers);
    $self->{_column_count} = scalar @$headers;
    if ($self->print_headers) {
        $self->output->print( join($self->separator, @$headers)."\n" );
    }

    if ( not defined $self->in_place_of_null_value ) {
        $self->in_place_of_null_value('');
    }

    return $self;
}

sub get_column_count {
    return $_[0]->{_column_count};
}

BEGIN {
    *Genome::Utility::IO::SeparatedValueWriter::print = \&Genome::Utility::IO::SeparatedValueWriter::write_one;
    *Genome::Utility::IO::SeparatedValueWriter::write = \&Genome::Utility::IO::SeparatedValueWriter::write_one;
}

sub write_one {
    my ($self, $data) = @_;

    $self->_validate_data_to_write($data)
        or return;

    return $self->output->print(
        join(
            $self->separator,
            map { defined $_ ? $_ : $self->in_place_of_null_value } map { $data->{$_} } @{$self->headers}
        )."\n"
    );
}

sub _validate_data_to_write {
    my ($self, $data) = @_;

    unless ( $data ) {
        $self->error_message("No data sent to 'write_one'");
        return;
    }

    my $reftype = Scalar::Util::reftype($data);
    unless ( $reftype and $reftype eq 'HASH' ) {
        $self->error_message("Need data as an hash ref to 'write_one'. Received:\n".Dumper($data));
        return;
    }

    unless ( %$data ) {
        $self->error_message("No data in data hash ref sent to 'write_one'");
        return;
    }

    my @headers = sort @{$self->headers};
    my @keys = sort keys %$data;
    unless ( @headers == @keys ) {
        # If we dont care about extra columns, all is well... unless we dont at least have the minimum required
        if ((!$self->ignore_extra_columns) || (@headers > @keys) ) {
                #  Bomb out if we dont want extra columns
                $self->error_message(
                    sprintf(
                        'Expected %d values, got %d in hash %s.',
                        scalar @headers,
                        scalar @keys,
                        Data::Dumper::Dumper($data),
                    )
                );
                return;
        }
    } else {
        unless ( Compare(\@headers,\@keys) ) {
            $self->error_message("Headers in data do not match headers being written:\n\tHEADERS:\n". Dumper(@headers) ."\tDATA_HASH_KEYS:\n". Dumper(@keys));
            return;
        }
    }

    return 1;
}

1;

=pod

=head1 Name

Genome::Utility::IO::SeparatedValueWriter

=head1 Synopsis

A stream based reader that splits each line by the given separator.  If no headers are given, they will be derived from the first line of the io, being split by the separator.

=head1 Usage

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

