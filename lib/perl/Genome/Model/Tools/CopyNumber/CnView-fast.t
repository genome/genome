#!/usr/bin/env genome-perl

use strict;
use warnings;
use above "Genome";
use Test::More tests => 5;

my $expected_results = $ENV{GENOME_TEST_INPUTS} . "/Genome-Model-Tools-CopyNumber-CnView-fast/2012-11-28";
ok(-d $expected_results, "test data dir is " . $expected_results)
  or die "cannot continue";

my $actual_results = Genome::Sys->create_temp_directory;
ok($actual_results, "test data writing to $actual_results")
  or die "cannot continue";

my $cmd = <<EOS;
  gmt copy-number cn-view \\
    --annotation-build=124434505 \\
    --cnv-file=/gscmnt/gc13001/info/model_data/2888915570/build129973671/variants/cnvs.hq \\
    --segments-file=/gscmnt/gc2013/info/model_data/2889110844/build130030495/PNC6/clonality/cnaseq.cnvhmm \\
    --output-dir=$actual_results \\
    --gene-targets-file=/gscmnt/sata132/techd/mgriffit/reference_annotations/GeneSymbolLists/CancerGeneCensusPlus_Sanger.txt \\
    --name='CancerGeneCensusPlus_Sanger' \\
    --chr=21 \\
EOS
eval { Genome::Sys->shellcmd(cmd => $cmd) };
ok(!$@, "no exceptions running the tool")
  or diag($@);

my @diff = `diff -r --brief -x '*.jpeg' $expected_results $actual_results`;
ok(@diff == 0, "no differences between expected and actual results") 
  or do {
      diag(@diff);
   };

my $file_check="-s $actual_results"."/CNView_CancerGeneCensusPlus_Sanger/chr21.jpeg";
ok($file_check, "chr21.jpeg found with non-zero size")
  or do {
      diag($file_check);
   };

if (@ARGV == 1 and $ARGV[0] eq 'KEEP') {
  my $stash = "/tmp/last-failed-cnview-fast-test";
  note("temp results moved to $stash");
  system "mv $actual_results $stash";
}

