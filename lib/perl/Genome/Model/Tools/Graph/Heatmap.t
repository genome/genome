#!/usr/bin/env genome-perl

use strict;
use warnings;

use above "Genome";
use Genome::Model::Tools::Graph::Heatmap;
use Test::More tests => 3;
use Digest::MD5 qw(md5_hex);
#plan "skip_all";

my $file = $ENV{GENOME_TEST_INPUTS} . "/Genome-Model-Tools-Graph/heatmap-test-matrix.csv";
my $tmp_dir = Genome::Sys->create_temp_directory('Genome-Model-Tools-Graph');
my $outfile = "$tmp_dir/heatmap-test-image.png";
my $columns = 3;
my $checksum = "007d3bb4cfa3bb2aacf152dcfa02aafa";

unlink $outfile;

my $hm = Genome::Model::Tools::Graph::Heatmap->create(
                                                         matrix => $file,
                                                         image => $outfile,
                                                         columns => $columns,
                                                        );
ok($hm->execute,'heatmap image generation');
ok((-e $outfile), 'output file exists');

#my $imagecontents = qx(cat $outfile);
#is(md5_hex($imagecontents), $checksum, 'image correct');
#ok(0,'firetest');

unlink $outfile;

my $nhfile = $ENV{GENOME_TEST_INPUTS} . "/Genome-Model-Tools-Graph/heatmap-noheader-matrix.csv";
my $hm_noheader = Genome::Model::Tools::Graph::Heatmap->create(
                                                         matrix => $nhfile,
                                                         image => $outfile,
                                                         columns => $columns,
                                                        );
my $ret;
eval { $ret = $hm_noheader->execute(); };
is($ret, undef, 'bad header');
unlink $outfile;

exit;

