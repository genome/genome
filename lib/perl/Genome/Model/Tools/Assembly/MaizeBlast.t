#!/usr/bin/env genome-perl

use strict;
use warnings;

use above "Genome";
use Test::More;

#TODO - find a sequence that'll finish blasting fast!

#for now just check to make sure maize bes blast db is there
foreach ('xnd', 'xni', 'xns', 'xnt') {
    ok (-s "/gscmnt/sata910/assembly/nthane/Maize_database/maize_454_db.$_", "Blast $_ file exists");
}

done_testing();
