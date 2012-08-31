package Genome::Model::Tools::Sx::IlluminaFastqReader;

use strict;
use warnings;

use Genome;

class Genome::Model::Tools::Sx::IlluminaFastqReader {
    is => 'Genome::Model::Tools::Sx::FastqReader',
};

sub read {
    my $seq = $_[0]->SUPER::read;
    return if not $seq;
    $seq->{qual} = join('', map { chr } map { ord($_) - 31 } split('', $seq->{qual}));
    return $seq;
}

1;

