package Genome::Model::Tools::Sx::PhredReader;

use strict;
use warnings;

use Genome;

class Genome::Model::Tools::Sx::PhredReader {
    is => 'Genome::Model::Tools::Sx::SeqReader',
    has => [ qual_file => { is => 'Text', is_optional => 1, }, ],
};

sub type { return 'phred'; }

sub create {
    my $class = shift;

    my $self = $class->SUPER::create(@_);
    return if not $self;

    $self->{_seq_reader} = Genome::Model::Tools::Sx::PhredSeqReader->create(
        file => $self->file,
    );
    return if not $self->{_seq_reader};

    return $self if not $self->qual_file;

    $self->{_qual_reader} = Genome::Model::Tools::Sx::PhredQualReader->create(
        file => $self->qual_file,
    );
    return if not $self->{_qual_reader};

    return $self;
}

sub read {
    my $self = shift;

    my $seq = $self->{_seq_reader}->read;
    return if not $seq;

    return $seq if not $self->{_qual_reader};

    my $qual = $self->{_qual_reader}->read;
    if ( not $qual ) {
        Carp::confess('No quality for sequence: '.$seq->{id});
    }
    if ( $seq->{id} ne $qual->{id} ) {
        Carp::confess('Fasta and quality ids do not match: '.$seq->{id}.' <=> '.$qual->{id});
    }
    if ( length($seq->{seq}) != length($qual->{qual}) ) {
        Carp::confess('Number of qualities does not match number of bases for fasta for '.$seq->{id}.'. Seq length: '.length($seq->{seq}).'. Qual length: '.length($qual->{qual}));
    }
    $seq->{qual} = $qual->{qual};

    return $seq;
}

sub close {
    my $self = shift;

    $self->{_seq_reader}->close if $self->{_seq_reader};
    $self->{_qual_reader}->close if $self->{_qual_reader};

    return 1;
}

1;

