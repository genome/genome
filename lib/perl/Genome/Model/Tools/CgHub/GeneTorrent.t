#! /gsc/bin/perl

use strict;
use warnings;

use above 'Genome';

use Test::More;

my $class = 'Genome::Model::Tools::CgHub::GeneTorrent';
use_ok($class) or die;

my $uuid = '387c3f70-46e9-4669-80e3-694d450f2919';
my $gene_torrent = Genome::Model::Tools::CgHub::GeneTorrent->create(
    uuid => $uuid,
    target_path => '/tmp',
);;
ok($gene_torrent, 'create gene torrent cmd');
my $rate_limit = $gene_torrent->rate_limit;
is($rate_limit, 10, 'correct rate limit');
is($gene_torrent->_build_command, "gtdownload --credential-file /gscuser/kochoa/mykey.pem --download https://cghub.ucsc.edu/cghub/data/analysis/download/$uuid --path /tmp --log stdout:verbose --verbose 2 --max-children 2 --rate-limit $rate_limit --inactivity-timeout ". (3 * 60 * 24), 'correct command');

done_testing();
