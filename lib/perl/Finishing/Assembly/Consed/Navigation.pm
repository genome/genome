package Finishing::Assembly::Consed::Navigation;

use strict;
use warnings;
no warnings 'reserved';

use Finfo::Std;

my %contig_name :name(contig_name:r) :isa(string);
my %desc :name(description:r) :isa(string);
my %read_name :name(read_name:o) :isa(string);
my %start :name(start:r) :isa('int non_neg');
my %stop :name(stop:r) :isa('int non_neg');
my %type :name(type:r) :isa('in_list CONSENSUS READ') :default('CONSENSUS');
my %acefile :name(acefile:o) :isa(string);

sub START
{
    my $self = shift;

    $self->error_msg("Navigation stop is less than the start")
        and return if $self->stop < $self->start;

    $self->error_msg("Type is READ, but no read name defined")
        and return if $self->type eq 'READ' and not defined $self->read_name;

    return 1;
}

sub length
{
    my $self = shift;

    return $self->stop - $self->start + 1;
}

1;

=pod

=head1 Name

Navigation

=head1 Synopsis

 Represents a navigation object for consed

=head1 Usage

use Navigation::Consed;

my $nav = Navigation::Consed->new
(
     contig_name => 'Contig44',
     description => 'Consensus Error',
     start => 455,
     stop => 455,
     type => 'CONSENSUS',
);

 or die Navigation::Consed->short_error_msg . "\n";

use Navigation::Consed::Reader;
use IO::File;

my $fh = IO::File->new("< list.txt")
    or die "$!\n";
my $reader = Navigation::Consed::Reader->new(io => $fh)
    or die;
while ( my $nav = $reader->next )
{
 ...
}

=head1 Accessors

=head2 contig_name

=head2 description

=head2 read_name

If type is 'READ', then this will be the name of the read.

=head2 start

=head2 stop

=head2 type

CONSENSUS or READ

=head1 Methods

=head2 length

The length between the start and stop

=head1 Author

Eddie Belter <ebelter@watson.wustl.edu>

=cut

#$HeadUR:$
#$Id: Navigation.pm 29586 2007-10-29 15:46:09Z ebelter $
