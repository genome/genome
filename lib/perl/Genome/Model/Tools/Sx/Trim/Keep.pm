package Genome::Model::Tools::Sx::Trim::Keep;

use strict;
use warnings;

use Genome;            

use Regexp::Common;

class Genome::Model::Tools::Sx::Trim::Keep {
    is => 'Genome::Model::Tools::Sx::Base',
    has => [
        length => {
            is => 'Integer',
            doc => 'The number of bases to keep',
        },    
     ],
};

sub help_brief {
    return 'Keep a set length of a sequence';
}

sub __errors__ {
    my $self = shift;
    my @errors = $self->SUPER::__errors__(@_);
    return @errors if @errors;
    if ( $self->length !~ /^$RE{num}{int}$/ or $self->length < 1 ) {
        push @errors, UR::Object::Tag->create(
            type => 'invalid',
            properties => [qw/ length /],
            desc => 'Keep length (length) is not a integer greater than 0 => '.$self->length,
        );
    }
    return @errors;
}

sub _eval_seqs {
    my ($self, $seqs) = @_;

    for my $seq (@$seqs) {
        $seq->{seq} = substr($seq->{seq}, 0, $self->length);
        $seq->{qual} = substr($seq->{qual},0, $self->length);
    }

    return $seqs;
}

1;

