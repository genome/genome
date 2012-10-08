package Genome::Model::Tools::Sx::Filter::Base;

use strict;
use warnings;

use Genome;

use Data::Dumper 'Dumper';

class Genome::Model::Tools::Sx::Filter::Base {
    is  => 'Genome::Model::Tools::Sx::Base',
    is_abstract => 1,
};

sub help_brief {
    return 'Filter sequences';
}

sub execute {
    my $self = shift;

    my $init = $self->_init;
    return if not $init;

    my $reader = $self->_reader;
    my $writer = $self->_writer;

    my $evaluator = $self->_create_sequence_evaluator;
    return if not $evaluator;

    SEQS: while ( my $seqs = $reader->read ) {
        next if not $evaluator->($seqs);
        $writer->write($seqs);
    }

    return 1;
}

sub _create_sequence_evaluator {
    my $self = shift;

    my @filters = $self->_create_filters;
    return if not @filters;

    return sub{
        for my $filter ( @filters ) {
            return if not $filter->($_[0]);
        }
        return 1;
    }
}

1;

