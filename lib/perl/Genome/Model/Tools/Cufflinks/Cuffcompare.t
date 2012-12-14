#!/usr/bin/env genome-perl

use strict;
use warnings;

use above 'Genome';

use File::Compare;
use Test::More tests => 10;

use_ok('Genome::Model::Tools::Cufflinks::Cuffcompare');

my $data_dir = "$ENV{GENOME_TEST_INPUTS}/Genome-Model-Tools-Cufflinks-Cuffcompare";

my $gtf_basename = 'test.gtf';
my $gtf_file = $data_dir .'/'. $gtf_basename;

my $fasta_file = $data_dir .'/test.fa';

my $expected_base = $data_dir .'/expected';
my $expected_combined_gtf_path = $expected_base .'.combined.gtf';
my $expected_loci_path = $expected_base .'.loci';
my $expected_stats_path = $expected_base .'.stats';
my $expected_tracking_path = $expected_base .'.tracking';

my $expected_refmap_path = $expected_base .'.'. $gtf_basename .'.refmap';
my $expected_tmap_path = $expected_base .'.'. $gtf_basename .'.tmap';

my $tmp_dir = File::Temp::tempdir('Cufflinks-Cuffcompare-'.Genome::Sys->username.'-XXXX',DIR => "$ENV{GENOME_TEST_TEMP}",CLEANUP => 1);
my $prefix = $tmp_dir .'/cuffcompare_test';

my $compare = Genome::Model::Tools::Cufflinks::Cuffcompare->create(
    use_version => '1.3.0',
    input_gtf_paths => $gtf_file,
    output_prefix => $prefix,
    reference_gtf_path => $gtf_file,
    reference_fasta_path => $fasta_file,
    include_contained => 1,
    generic_gtf_input => 1,
);
isa_ok($compare,'Genome::Model::Tools::Cufflinks::Cuffcompare');
ok($compare->execute,'execute command '. $compare->command_name);

#This file appears to be different when run as another user or maybe on a different machine... need to figure this out
#ok( (compare($compare->combined_gtf_path,$expected_combined_gtf_path) == 0),'combined transcripts are identical');
ok( -s $compare->combined_gtf_path, 'combined gtf file exists');

ok( (compare($compare->loci_path,$expected_loci_path) == 0),'loci are identical');
#The path to the files is included in the output file... needs scrubbed or something
#ok( (compare($compare->stats_path,$expected_stats_path) == 0),'stats file'. $compare->stats_path .' is identical to expected '. $expected_stats_path);
ok( (compare($compare->tracking_path,$expected_tracking_path) == 0),'tracking are identical');

my @refmap_paths = @{$compare->refmap_paths};
is(scalar(@refmap_paths),1,'Found 1 refmap path');
ok( (compare($refmap_paths[0],$expected_refmap_path) == 0),'refmap are identical');

my @tmap_paths = @{$compare->tmap_paths};
is(scalar(@tmap_paths),1,'Found 1 refmap path');
ok( (compare($tmap_paths[0],$expected_tmap_path) == 0),'tmap are identical');
