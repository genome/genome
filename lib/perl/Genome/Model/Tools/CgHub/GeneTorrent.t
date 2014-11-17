#! /gsc/bin/perl

use strict;
use warnings;

use above 'Genome';

use Test::More;

use_ok('Genome::Model::Tools::GeneTorrent') or die;

my $gene_torrent = Genome::Model::Tools::GeneTorrent->create;
ok($gene_torrent, 'create gene torrent cmd');
is($gene_torrent->rate_limit, 10, 'correct rate limit');

done_testing();
