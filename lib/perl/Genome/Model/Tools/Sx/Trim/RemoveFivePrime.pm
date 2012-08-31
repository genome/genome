package Genome::Model::Tools::Sx::Trim::RemoveFivePrime;

use strict;
use warnings;

use Genome;            

use Regexp::Common;

class Genome::Model::Tools::Sx::Trim::RemoveFivePrime {
    is => 'Genome::Model::Tools::Sx::Base',
    has => [
        length => {
            is => 'Integer',
            doc => 'The number of bases to remove.',
        },    
     ],
};

sub help_brief {
    return 'Remove bases from the 5 prime (left) end';
}

sub __errors__ {
    my $self = shift;
    my @errors = $self->SUPER::__errors__(@_);
    return @errors if @errors;
    if ( $self->length !~ /^$RE{num}{int}$/ or $self->length < 1 ) {
        push @errors, UR::Object::Tag->create(
            type => 'invalid',
            properties => [qw/ length /],
            desc => 'The remove length (length) is not a integer greater than 0 => '.$self->length,
        );
    }
    return @errors;
}

sub _eval_seqs {
    my ($self, $seqs) = @_;

    for my $seq ( @$seqs ) {
        if ( $self->length > length($seq->{seq}) ) {
            $seq->{seq} = '';
            $seq->{qual} = '';
        }
        else {
            $seq->{seq} = substr($seq->{seq}, $self->length);
            $seq->{qual} =substr($seq->{qual}, $self->length);
        }
    }

    return $seqs;
}

1;

