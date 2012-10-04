package Genome::Model::Tools::Sx::Filter::ByMinNonNs;

use strict;
use warnings;

use Genome;            

use Regexp::Common;

class Genome::Model::Tools::Sx::Filter::ByMinNonNs {
    is => 'Genome::Model::Tools::Sx::Filter::Base',
    has => [
        minimum => {
            is => 'Number',
            doc => 'The minimum number of non Ns in a sequence.',
        }    
    ],
};

sub help_synopsis {
    return <<HELP
    Filter sequences by minimum non Ns. If one sequence in a set fails, the entire set is discarded.
HELP
}

sub help_brief {
    return 'Filter sequences by minimum non Ns';
}

sub __errors__ {
    my $self = shift;

    my @errors = $self->SUPER::__errors__(@_);
    return @errors if @errors;

    my $minimum = $self->minimum;
    if ( defined $minimum and ( $minimum !~ /^$RE{num}{int}$/ or $minimum < 0 ) ) {
        push @errors, UR::Object::Tag->create(
            type => 'invalid',
            properties => [qw/ minimum /],
            desc => "Minimum number of non Ns ($minimum) must be a greater than or equal to 0",
        );
    }

    return @errors;
}

sub _create_filters {
    my $self = shift;

    my $minimum = $self->minimum;
    return sub{
        my $seqs = shift;
        for my $seq ( @$seqs ) {
            my $ns = () = $seq->{seq} =~ /n/ig;
            if ( (length($seq->{seq}) - $ns) < $minimum ) {
                return;
            }
        }
        return 1;
    }
}

1;

