package Genome::Model::Tools::Sx::PhredSeqReader;

use strict;
use warnings;

use Genome;

class Genome::Model::Tools::Sx::PhredSeqReader {
    is => 'Genome::Model::Tools::Sx::PhredReaderBase',
};

sub read {
    my $self = shift;

    my %seq;
    @seq{qw/ id desc seq /} = $self->_parse_io($self->{_file});
    return if not $seq{seq};
    $seq{seq} =~ tr/ \t\n\r//d;	# Remove whitespace

    return \%seq;
}

1;

