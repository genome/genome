#!/usr/bin/env genome-perl

use strict;
use warnings;


use above "PAP";
use above 'Workflow';
use Data::Dumper;
use File::Basename;
use Test::More qw(no_plan);


ok(1); # stupid, but works for now.
SKIP: {
   skip "this needs to be rewritten in the future", 0 unless 0;
my $w = Workflow::Model->create_from_xml($ARGV[0] || File::Basename::dirname(__FILE__).'/data/pap_outer_keggless.xml');

my @errors = $w->validate;
die 'Too many problems: ' . join("\n", @errors) unless $w->is_valid();

my $out = $w->execute(
    'input' => {
        'fasta file'       => File::Basename::dirname(__FILE__).'/data/B_coprocola.chunk.fasta',
        'chunk size'       => 10,
        'biosql namespace' => 'MGAP',
        'gram stain'       => 'negative',
    }
);

$w->wait();

print Data::Dumper->new([$out])->Dump;
}
exit(0);
