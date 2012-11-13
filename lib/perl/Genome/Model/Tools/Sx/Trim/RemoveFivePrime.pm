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
    my $length = $self->length;
    if ( $length !~ /^$RE{num}{int}$/ or $length < 1 ) {
        push @errors, UR::Object::Tag->create(
            type => 'invalid',
            properties => [qw/ length /],
            desc => 'The remove length (length) is not a integer greater than 0 => '.$length,
        );
    }
    return @errors;
}

sub _create_evaluator {
    my $self = shift;

    my $length = $self->length;
    return sub{
        for my $seq ( @{$_[0]} ) {
            if ( $length > length($seq->{seq}) ) {
                $seq->{seq} = '';
                $seq->{qual} = '';
            }
            else {
                $seq->{seq} = substr($seq->{seq}, $length);
                $seq->{qual} =substr($seq->{qual}, $length);
            }
        }
        return 1;
    }
}

1;

