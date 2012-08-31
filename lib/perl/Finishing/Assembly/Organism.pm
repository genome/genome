package Finishing::Assembly::Organism;

use strict;
use warnings;

use base 'Finishing::Assembly::Item';

#- NAMES -#
sub spaced_name { return species_name(@_); }
sub species_name
{
    my $self = shift;

    my @tokens = split(/\_/, $self->name);
    $tokens[0] = ucfirst $tokens[0];

    return join(' ', ucfirst($tokens[0]), @tokens[1..$#tokens]);
}

1;

#$HeadURL$
#$Id$
