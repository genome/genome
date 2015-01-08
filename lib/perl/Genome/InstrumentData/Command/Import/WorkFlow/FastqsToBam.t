#!/usr/bin/env genome-perl

BEGIN {
    $ENV{UR_DBI_NO_COMMIT} = 1;
    $ENV{UR_COMMAND_DUMP_STATUS_MESSAGES} = 1;
    $ENV{UR_COMMAND_DUMP_STATUS_MESSAGES} = 1;
}

use strict;
use warnings;

use above 'Genome';

require Genome::Utility::Test;
require File::Compare;
require File::Spec;
require File::Temp;
use Test::More;

use_ok('Genome::InstrumentData::Command::Import::WorkFlow::FastqsToBam') or die;
my $test_dir = Genome::Utility::Test->data_dir_ok('Genome::InstrumentData::Command::Import', 'fastq/v2') or die;
use_ok('Genome::InstrumentData::Command::Import::WorkFlow::Helpers') or die;
my $helpers = Genome::InstrumentData::Command::Import::WorkFlow::Helpers->get;

my $tmp_dir = File::Temp::tempdir(CLEANUP => 1);

my @source_fastq_base_names = (qw/ input.1.fastq.gz input.2.fastq /);
my @source_fastq_paths;
for my $source_fastq_base_name ( @source_fastq_base_names ) {
    my $source_fastq_path = $tmp_dir.'/'.$source_fastq_base_name;
    push @source_fastq_paths, $source_fastq_path;
    Genome::Sys->create_symlink($test_dir.'/'.$source_fastq_base_name, $source_fastq_path);
}

my $sample_name = '__TEST_SAMPLE__';
my $library_name = join('-', $sample_name, 'extlibs');

my $cmd = Genome::InstrumentData::Command::Import::WorkFlow::FastqsToBam->execute(
    working_directory => $tmp_dir,
    fastq_paths => \@source_fastq_paths,
    sample_name => $sample_name,
    library_name => $library_name,
);
ok($cmd->result, 'execute');
my $output_bam_path = $cmd->output_bam_path;
is($output_bam_path, $tmp_dir.'/__TEST_SAMPLE__.bam', 'bam path named correctly');
ok(-s $output_bam_path, 'bam path exists');
my $expected_bam = File::Spec->join($test_dir, 'fastqs-to-bam.bam');
is(File::Compare::compare($output_bam_path, $expected_bam), 0, 'bam matches');
is(File::Compare::compare($output_bam_path.'.flagstat', $expected_bam.'.flagstat'), 0, 'flagstat matches');

for ( my $i = 0; $i < @source_fastq_paths; $i++ ) {
    ok(!glob($source_fastq_paths[$i].'*'), 'removed fastq '.($i+1).' after conversion to bam');
}

#print "$tmp_dir\n"; <STDIN>;
done_testing();
