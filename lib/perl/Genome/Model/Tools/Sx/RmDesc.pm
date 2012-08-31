package Genome::Model::Tools::Sx::RmDesc;

use strict;
use warnings;

use Genome;

class Genome::Model::Tools::Sx::RmDesc {
    is  => 'Genome::Model::Tools::Sx::Base',
};

sub help_brief {
    return 'Remove the description from sequences';
}

sub execute {
    my $self = shift;

    my $init = $self->_init;
    return if not $init;

    my $reader = $self->_reader;
    my $writer = $self->_writer;

    while ( my $seqs = $reader->read ) {
        for my $seq ( @$seqs ) { $seq->{desc} = undef }
        $writer->write($seqs);
    }

    return 1;
}

1;

