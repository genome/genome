package Genome::Model::Tools::Sx::IlluminaFastqWriter;

use strict;
use warnings;

use Genome;

class Genome::Model::Tools::Sx::IlluminaFastqWriter {
    is => 'Genome::Model::Tools::Sx::FastqWriter',
};

sub write {
    my ($self, $seq) = @_;
    # sanger to illumina
    $seq->{qual} = join('', map { chr } map { ord($_) + 31 } split('', $seq->{qual}));
    return $self->SUPER::write($seq);
}

1;

