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

my $tmpdir = File::Temp::tempdir('GMT-Pindel-RunPindel-XXXXX', CLEANUP => 1, TMPDIR => 1);
my $test_data = Genome::Config::get('test_inputs') . "/Genome-Model-Tools-Pindel-RunPindel";
my $tumor_bam = $test_data."/true_positive_tumor_validation.bam";
my $normal_bam = $test_data."/true_positive_normal_validation.bam";
my $expected_directory = $test_data ."/expected_2";
my $expected_indels_hq = $expected_directory."/10/indels.hq";
my $actual_indels_hq = $tmpdir."/10/indels.hq";

my $pindel = Genome::Model::Tools::Pindel::RunPindel->create(
                aligned_reads_input => $tumor_bam,
                control_aligned_reads_input => $normal_bam,
                reference_build_id => $refbuild_id,
                version => '0.5',
                output_directory => $tmpdir,
                chromosome => '10', );


ok($pindel, 'run-pindel command created');

my $result = $pindel->execute;
is($result, 1, 'Testing for execution.  Expecting 1.  Got: '.$result);

is(compare($actual_indels_hq,$expected_indels_hq),0,'Output for v0.5 is identical to expected output');

done_testing();
