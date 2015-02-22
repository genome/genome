#!/usr/bin/env genome-perl

use strict;
use warnings;

use above 'Genome';

require Genome::Utility::Test;
use Test::Exception;
use Test::More;

use_ok('Genome::Model::Tools::CgHub::GeneTorrent') or die;

my $data_dir = Genome::Utility::Test->data_dir_ok('Genome::Model::Tools::CgHub');
my $path = File::Spec->join($data_dir, 'exists.bam');

# So we don't actually send a request
my $run_command_cnt = 0;
sub Genome::Model::Tools::CgHub::GeneTorrent::_run_command { $run_command_cnt++; return 1; };

# Success
my $uuid = '387c3f70-46e9-4669-80e3-694d450f2919';
my $gene_torrent = Genome::Model::Tools::CgHub::GeneTorrent->create(uuid => $uuid, target_path => $path);
ok($gene_torrent, 'create gene torrent cmd');
is($gene_torrent->source_url, 'https://cghub.ucsc.edu/cghub/data/analysis/download/'.$uuid, 'correct source_url');
is($gene_torrent->rate_limit, 10, 'correct rate limit');
is(
    $gene_torrent->_build_command,
    'gtdownload'
    . ' --credential-file '.$gene_torrent->credential_file
    . ' --download '. $gene_torrent->source_url
    . ' --path ' . $gene_torrent->target_path
    . ' --log stdout:verbose'
    . ' --verbose 2'
    . ' --max-children 2'
    . ' --rate-limit '.$gene_torrent->rate_limit # mega-BYTES per second (see internet_download_mbps above)
    . ' --inactivity-timeout ' . 3 * 60 * 24   # in minutes - instead of bsub -W
    ,
    'correct cmd',
);
ok($gene_torrent->execute, 'execute ok');

# Fail
$gene_torrent = Genome::Model::Tools::CgHub::GeneTorrent->create(uuid => $uuid, target_path => 'does_not_exist');
ok($gene_torrent, 'create gene torrent cmd');
throws_ok(sub{ $gene_torrent->execute; }, qr//, 'execute fails b/c target path does not exist');

is($run_command_cnt, 2, 'run command was invoked correctly');

done_testing();
