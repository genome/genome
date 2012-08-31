#!/usr/bin/env genome-perl

use strict;
use warnings;

use Genome;

main();
exit;

sub main {
    my @sr = GSC::Equipment::Solexa::Run->get( run_status => 'active' );
    Genome::Search->add(@sr);
}
