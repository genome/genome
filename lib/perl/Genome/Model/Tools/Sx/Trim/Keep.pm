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
    my $length = $self->length;
    if ( $length !~ /^$RE{num}{int}$/ or $length < 1 ) {
        push @errors, UR::Object::Tag->create(
            type => 'invalid',
            properties => [qw/ length /],
            desc => 'Keep length (length) is not a integer greater than 0 => '.$length,
        );
    }
    return @errors;
}

sub _create_evaluator {
    my $self = shift;

    my $length = $self->length;
    return sub{
        for my $seq ( @{$_[0]} ) {
            $seq->{seq} = substr($seq->{seq}, 0, $length);
            $seq->{qual} = substr($seq->{qual},0, $length);
        }
        return 1;
    }
}

1;

