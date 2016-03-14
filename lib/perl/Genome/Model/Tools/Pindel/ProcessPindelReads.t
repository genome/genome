#!/usr/bin/env genome-perl

use strict;
use warnings;

use File::Path;
use File::Temp;
use File::Compare;
use Test::More;
use above 'Genome';

my $archos = `uname -a`;
if ($archos !~ /64/) {
    plan skip_all => "Must run from 64-bit machine";
}

# Caching refseq in /var/cache/tgi-san. We gotta link these files to a tmp dir for tests so they don't get copied
my $refbuild_id = 101947881;
my $ref_seq_build = Genome::Model::Build::ImportedReferenceSequence->get($refbuild_id);
ok($ref_seq_build, 'human36 reference sequence build') or die;
my $refseq_tmp_dir = File::Temp::tempdir(CLEANUP => 1);
no warnings;
*Genome::Model::Build::ReferenceSequence::local_cache_basedir = sub { return $refseq_tmp_dir; };
*Genome::Model::Build::ReferenceSequence::copy_file = sub { 
    my ($build, $file, $dest) = @_;
    symlink($file, $dest);
    is(-s $file, -s $dest, 'linked '.$dest) or die;
    return 1; 
};
use warnings;

my $tmpdir = File::Temp::tempdir('GMT-Pindel-ProcessPindelReads-XXXXX', CLEANUP => 1, TMPDIR => 1);
my $test_data = Genome::Config::get('test_inputs') . "/Genome-Model-Tools-Pindel-ProcessPindelReads";
my $input_02 = "$test_data/indels_all_sequences.0.2";
my $input_04 = "$test_data/indels_all_sequences.0.4";
my $input_05 = "$test_data/indels_all_sequences.0.5";
my $output_02 = "$tmpdir/indels.hq.v02.bed";
my $output_04 = "$tmpdir/indels.hq.v04.bed";
my $output_05 = "$tmpdir/indels.hq.v05.bed";
my $output_06 = "$tmpdir/indels.hq.v06.bed";
my $big_output_05 = "$tmpdir/indels.hq.v05.bed.big_insertions";
my $expected_output_02 = "$test_data/expected/indels.hq.v02.bed";
my $expected_output_04 = "$test_data/expected/indels.hq.v04.bed";
my $expected_output_05 = "$test_data/expected/indels.hq.v05.bed";
my $expected_output_06 = "$test_data/expected/indels.hq.v06.bed";

# Test Pindel v0.5 output
my $ppr_cmd_06 = Genome::Model::Tools::Pindel::ProcessPindelReads->create(
                input_file => $input_05,
                output_file => $output_06,
                big_output_file => $big_output_05,
                reference_build_id => $refbuild_id,
                mode => 'to_bed', );

ok($ppr_cmd_06, 'process-pindel-reads command created');

my $result = $ppr_cmd_06->execute;
is($result, 1, 'Testing for execution.  Expecting 1.  Got: '.$result);

is(compare($output_06,$expected_output_06),0,'Output for v0.5 is identical to expected output');

done_testing();
