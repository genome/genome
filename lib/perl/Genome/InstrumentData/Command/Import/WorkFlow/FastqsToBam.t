#! /gsc/bin/perl

use strict;
use warnings;

use above 'Genome';

require Genome::Utility::Test;
require File::Compare;
require File::Temp;
use Test::More;

use_ok('Genome::InstrumentData::Command::Import::WorkFlow::VerifyMd5') or die;
my $test_dir = Genome::Utility::Test->data_dir_ok('Genome::InstrumentData::Command::Import', 'v1') or die;
my $tmp_dir = File::Temp::tempdir(CLEANUP => 1);

my @source_fastq_base_names = (qw/ input.1.fastq.gz input.2.fastq /);
my @source_fastq_paths = map { $tmp_dir.'/'.$_ } @source_fastq_base_names;
for my $source_fastq_base_name ( @source_fastq_base_names ) {
    my $source_fastq_path = $tmp_dir.'/'.$source_fastq_base_name;
    push @source_fastq_paths, $source_fastq_path;
    Genome::Sys->create_symlink($test_dir.'/'.$source_fastq_base_name, $source_fastq_path);
}

my $sample = Genome::Sample->__define__(id => -1, name => '__TEST_SAMPLE__');
ok($sample, 'define sample');

my $cmd = Genome::InstrumentData::Command::Import::WorkFlow::FastqsToBam->execute(
    working_directory => $tmp_dir,
    fastq_paths => \@source_fastq_paths,
    sample => $sample,
);
ok($cmd->result, 'execute');
my $bam_path = $cmd->bam_path;
is($bam_path, $tmp_dir.'/__TEST_SAMPLE__.bam', 'bam path named correctly');
ok(-s $bam_path, 'bam path exists');
is(File::Compare::compare($bam_path, $test_dir.'/input.fastq.unsorted.bam'), 0, 'bam matches');

#print "$tmp_dir\n"; <STDIN>;
done_testing();
