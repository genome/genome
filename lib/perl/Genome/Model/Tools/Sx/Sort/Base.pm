package Genome::Model::Tools::Sx::Sort::Base;

use strict;
use warnings;

use Genome;

require IPC::Open2;

class Genome::Model::Tools::Sx::Sort::Base {
    is  => 'Genome::Model::Tools::Sx::Base',
};

sub help_brief {
    return 'Sort sequences by average quality';
}

sub execute {
    my $self = shift;

    my $init = $self->_init;
    return if not $init;

    my $pid = IPC::Open2::open2(\*SORTED, \*UNSORTED, 'sort', $self->_sort_params);
    if ( not $pid ) {
        $self->error_message('Failed to create the sort command via open2');
        return;
    }

    my $flatten = $self->can('_flatten');
    my $inflate = $self->can('_inflate');

    my $reader = $self->_reader;
    while ( my $seqs = $reader->read ) {
        print UNSORTED ($flatten->(@$seqs), "\n");
    }
    close UNSORTED;

    my $writer = $self->_writer;
    while ( my $line = <SORTED> ) {
        $writer->write( $inflate->($line) );
    }
    close SORTED;
    waitpid($pid, 0);

    return 1;
}

1;

