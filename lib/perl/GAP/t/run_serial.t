#!/usr/bin/env genome-perl

use strict;
use warnings;

use above 'Workflow';
use above 'GAP';
use Data::Dumper;
use File::Basename;
use Test::More qw(no_plan);

my $w = Workflow::Model->create_from_xml(File::Basename::dirname(__FILE__).'/data/repeatmasker_outer.xml');

ok($w);
print join("\n", $w->validate) . "\n";

my $out = $w->execute(
    'input' => {
        'fasta_file' => File::Basename::dirname(__FILE__).'/data/C_elegans.ws184.fasta.bz2',
    }
);

$w->wait;

print Data::Dumper->new([$out])->Dump;

