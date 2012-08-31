#!/usr/bin/env genome-perl

use strict;
use warnings;
use above 'Genome';
use Test::More tests => 2;


use_ok( 'Genome::Model::Tools::ListReferenceSequences' );


my $result = Genome::Model::Tools::ListReferenceSequences->execute;

ok($result, 'Test ran ok');


