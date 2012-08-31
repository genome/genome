package Genome::Model::Tools::Sx::Filter::ByMinLength;

use strict;
use warnings;

use Genome;            

use Regexp::Common;

class Genome::Model::Tools::Sx::Filter::ByMinLength {
    is => 'Genome::Model::Tools::Sx::Filter::Base',
    has => [
        length => {
            is => 'Number',
            doc => 'The minimum length of the sequences.',
        }    
    ],
};

sub help_synopsis {
    return <<HELP
    Filter fastq sequences by length. Considers all sequences in the set.
HELP
}

sub help_brief {
    return 'Filter sequences by minimum length';
}

sub __errors__ {
    my $self = shift;

    my @errors = $self->SUPER::__errors__(@_);
    return @errors if @errors;

    my $length = $self->length;
    if ( not defined $length or $length !~ /^$RE{num}{int}$/ or $length < 1 ) {
        push @errors, UR::Object::Tag->create(
            type => 'invalid',
            properties => [qw/ length /],
            desc => "Length ($length) must be a positive integer greater than 1",
        );
    }

    return @errors;
}

sub _create_filters {
    my $self = shift;

    my $length = $self->length;
    return sub{
        my $seqs = shift;
        for my $seq ( @$seqs ) {
            unless ( length $seq->{seq} >= $length ) {
                return;
            }
        }
        return 1;
    }
}

1;

