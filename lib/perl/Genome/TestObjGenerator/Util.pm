package Genome::TestObjGenerator::Util;

use strict;
use warnings;
use Genome;

my %count;
sub generate_name {
    my $base_name = shift;
    $count{$base_name}++;
    return join('_', $base_name, $count{$base_name});
}

1;
