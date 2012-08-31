package Genome::Model::Tools::Sx::PhredEnhancedSeqReader;

# TODO
# Add functionality to put in the IUB code for ambiguous bases

use strict;
use warnings;

use Genome;

class Genome::Model::Tools::Sx::PhredEnhancedSeqReader {
    is => 'Genome::Model::Tools::Sx::PhredReaderBase',
};

sub read {
    my $self = shift;

    my %seq;
    @seq{qw/ id desc seq /} = $self->_parse_io($self->{_file});
    return if not $seq{seq};
    $seq{seq} =~ tr/ \t\n\r//d;	# Remove whitespace
    $seq{seq} =~ s#{(\w*),([\w,])*}#$1#g; # Use the first base

    return \%seq;
}

1;

