#:boberkfe is this used anymore?

package Genome::Model::AlignedBase;

use strict;
use warnings;

use constant alphabet_indexes => {
                                    '-' => 0,
                                    'A' => 1,
                                    'C' => 2,
                                    'G' => 3,
                                    'T' => 4,
                                };

sub new{
    my ($pkg, %params) = @_;
    
    my $self = {
        alignment_probability   => undef,
        base_probability_vector => undef,
    };

    $self = { %$self, %params  };
    
    return bless $self, $pkg;
}



1;
