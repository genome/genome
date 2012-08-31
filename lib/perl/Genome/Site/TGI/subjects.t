#!/usr/bin/env genome-perl

use strict;
use warnings;
use above "Genome";
use Test::More;

my @classes = qw/
    DNAResource           
    DNAResourceItem       
    IPRProduct            
    SolexaRun
    GenomicDNA            
    Ligation              
    PCRProduct            
/;

plan tests => scalar(@classes) * 3;

for my $class (@classes) {
    my $c2 = 'Genome::Site::TGI::'.$class;
    my $o = eval {
        my $i = $c2->create_iterator();
        my $o = $i->next;
        return $o;
    };
    ok($o, "got object for $c2");
    ok($o->id, "got id " . $o->id);
    ok($o->name, "got name " . $o->name);
}

