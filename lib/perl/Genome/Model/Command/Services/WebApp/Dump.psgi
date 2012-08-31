use strict;
use warnings;

use Data::Dumper;

sub {
    my ($env) = @_;
    return [
        200,
        [ 'Content-type' => 'text/plain' ],
        [ Data::Dumper->new( [$env] )->Dump ]
    ];
};
