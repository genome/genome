package Genome::TestObjGenerator::Util;

use strict;
use warnings;
use Genome;

my $name_count = 0;

sub generate_name {
    my $base_name = shift;
    $name_count++;
    return $base_name."_".$name_count; 
}

1;

