package Genome::Model::Tools::Sx::Filter::ByMaxNs;

use strict;
use warnings;

use Genome;            

use Regexp::Common;

class Genome::Model::Tools::Sx::Filter::ByMaxNs {
    is => 'Genome::Model::Tools::Sx::Filter::Base',
    has => [
        maximum => {
            is => 'Number',
            doc => 'The maximum number of Ns allowed in a sequence.',
        }    
    ],
};

sub help_synopsis {
    return <<HELP
    Filter sequences by maximum number of Ns in a sequence. If one sequence in a set fails, the entire set is discarded.
HELP
}

sub help_brief {
    return 'Filter sequences by maximum number of Ns';
}

sub __errors__ {
    my $self = shift;

    my @errors = $self->SUPER::__errors__(@_);
    return @errors if @errors;

    my $maximum = $self->maximum;
    if ( defined $maximum and ( $maximum !~ /^$RE{num}{int}$/ or $maximum < 0 ) ) {
        push @errors, UR::Object::Tag->create(
            type => 'invalid',
            properties => [qw/ maximum /],
            desc => "Maximum number of Ns ($maximum) must be a greater than or equal to 0",
        );
    }

    return @errors;
}

sub _create_filters {
    my $self = shift;

    my $maximum = $self->maximum;
    return sub{
        my $seqs = shift;
        for my $seq ( @$seqs ) {
            my $ns = () = $seq->{seq} =~ /n/ig;
            if ( $ns > $maximum ) {
                return;
            }
        }
        return 1;
    }
}

1;

