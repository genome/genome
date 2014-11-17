#! /gsc/bin/perl

use strict;
use warnings;

use above 'Genome';

use Test::More;

use_ok('Genome::Model::Tools::CgHub::GeneTorrent') or die;

my $gene_torrent = Genome::Model::Tools::CgHub::GeneTorrent->create(uuid => 'BLAH', target_path => 'PATH');
ok($gene_torrent, 'create gene torrent cmd');
is($gene_torrent->rate_limit, 10, 'correct rate limit');
is(
    $gene_torrent->_build_command,
    'gtdownload'
    . ' --credential-file /gscuser/kochoa/mykey.pem'    # TODO: do not hardcode
    . ' --download https://cghub.ucsc.edu/cghub/data/analysis/download/' . $gene_torrent->uuid
    . ' --path ' . $gene_torrent->target_path
    . ' --log stdout:verbose'
    . ' --verbose 2'
    . ' --max-children 2'
    . ' --rate-limit '.$gene_torrent->rate_limit # mega-BYTES per second (see internet_download_mbps above)
    . ' --inactivity-timeout ' . 3 * 60 * 24   # in minutes - instead of bsub -W
    ,
    'correct cmd',
);

done_testing();
