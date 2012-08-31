package Genome::Model::Tools::MetagenomicClassifier::ClassificationReader;

use strict;
use warnings;

use Genome;

class Genome::Model::Tools::MetagenomicClassifier::ClassificationReader {
    has => [
        file => {
            is => 'Text',
            doc => 'File of classifications to read.',
        },
        _io => { },
    ],
};

sub create {
    my $class = shift;

    my $self = $class->SUPER::create(@_);
    return if not $self;

    if ( not $self->file ) {
        $self->error_message('No file given to classification reader');
        return;
    }

    my $fh = eval { Genome::Sys->open_file_for_reading($self->file) };
    if ( not $fh ) {
        $self->error_message('Failed to open file ('.$self->file.') for reading: '.$@);
        return;
    }
    $self->_io($fh);

    return $self;
}

sub read {
    my $self = shift;

    my $line = $self->_io->getline;
    return if not $line;

    chomp $line;
    my ($id, $complemented, @taxa_and_confs) = split(';', $line);
    if ( @taxa_and_confs < 6 ) {
        Carp::confess('Failed to parse line in hmp_fix_ranks. Expected an id, orientation and at least 6 taxa: '.$line);
    }

    my %classification = (
        id => $id,
        complemented => ( $complemented eq '-' ? 1 : 0 ),
    );

    my @ranks = (qw/ root domain phylum class order family genus species /);
    for ( my $i = 0; $i <= $#taxa_and_confs; $i++ ) {
        my %taxon;
        $classification{$ranks[$i]} = \%taxon;
        next if not $taxa_and_confs[$i];
        @taxon{qw/ id confidence /} = split(':', $taxa_and_confs[$i]);
    }

    return \%classification;
}

1;

