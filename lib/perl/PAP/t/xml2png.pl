#!/usr/bin/env genome-perl

use strict;
use warnings;


use lib '/gscmnt/temp212/bioperl-svn/bioperl-live';
use lib '/gscmnt/temp212/bioperl-svn/bioperl-run';

use above 'Workflow';

use Data::Dumper;

my $w = Workflow::Model->create_from_xml($ARGV[0] || 'data/pap_outer.xml');

print join("\n", $w->validate) . "\n";

print $w->as_png($ARGV[1] || 'pap.png');
