#!/usr/bin/env genome-perl
use strict;
use warnings;

use Finishing::Assembly::Factory; 
my $fac = Finishing::Assembly::Factory->connect("cmap_user");
my $assembly = $fac->get_organism("pan_troglodytes")->get_assembly("2.1_051011");
my $counter = 0;
for( 1..1050){
    print ++$counter;
    my $contig = $assembly->get_contig("Contig9.$_");
    print $contig->name."\n";
}
