#!/usr/bin/env genome-perl

use strict;
use warnings;
use above "Genome";
use Test::More tests => 6;
use Genome::Model::Tools::CopyNumber::CnView;
use Genome::Utility::Test;
use File::Spec;

my $expected_results = Genome::Config::get('test_inputs') . "/Genome-Model-Tools-CopyNumber-CnView-slow/2012-12-03";
ok(-d $expected_results, "test data dir is " . $expected_results)
  or die "cannot continue";

my $actual_results = Genome::Sys->create_temp_directory;
ok($actual_results, "test data writing to $actual_results")
  or die "cannot continue";

my $test_data_dir = Genome::Utility::Test->data_dir('Genome::Model::Tools::CopyNumber::CnView', '2016-05-10_test-data');

my $cnv_file = File::Spec->join($test_data_dir, 'cnvs.hq');
my $segments_file = File::Spec->join($test_data_dir, 'cnaseq.cnvhmm');

my $db_name = 'tgi/cancer-annotation/human/build37-20130401.1';
my $db = Genome::Db->get($db_name);
my $targets_file = $db->data_set_path('GeneSymbolLists/CancerGeneCensusPlus_Sanger.txt');

my $cmd = <<EOS;
  gmt copy-number cn-view \\
    --annotation-build=124434505 \\
    --cnv-file=$cnv_file \\
    --segments-file=$segments_file \\
    --output-dir=$actual_results \\
    --gene-targets-file=$targets_file \\
    --name='CancerGeneCensusPlus_Sanger' \\
    --verbose \\
    --cancer-annotation-db $db_name \\
EOS
eval { Genome::Sys->shellcmd(cmd => $cmd) };
ok(!$@, "no exceptions running the tool")
  or diag($@);

my @diff = `diff -r --brief -x '*.jpeg' -x '*.stdout' -x '*.stderr' $expected_results $actual_results`;
ok(@diff == 0, "no differences between expected and actual results") 
  or do {
      diag(@diff);
   };

my $file_check1 = "$actual_results"."/CNView_CancerGeneCensusPlus_Sanger/Gains_chr21.jpeg";
ok(-s $file_check1, "Gains_chr21.jpeg found with non-zero size")
  or do {
      diag($file_check1);
   };

my $file_check2 = "$actual_results"."/CNView_CancerGeneCensusPlus_Sanger/Losses_chr21.jpeg";
ok(-s $file_check2, "Losses_chr21.jpeg found with non-zero size")
  or do {
      diag($file_check2);
   };

if (@ARGV == 1 and $ARGV[0] eq 'KEEP') {
  my $stash = "/tmp/last-failed-cnview-slow-test";
  note("temp results moved to $stash");
  if (-e $stash){
    system("rm -fr $stash");
  }
  system "mv $actual_results $stash";
}
