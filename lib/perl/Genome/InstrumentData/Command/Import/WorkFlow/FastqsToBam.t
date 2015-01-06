#!/usr/bin/env genome-perl

use strict;
use warnings;

use above 'Genome';

require Genome::Utility::Test;
require File::Compare;
require File::Temp;
use Test::More;

use_ok('Genome::InstrumentData::Command::Import::WorkFlow::FastqsToBam') or die;
my $test_dir = Genome::Utility::Test->data_dir_ok('Genome::InstrumentData::Command::Import', 'fastq/v1') or die;
use_ok('Genome::InstrumentData::Command::Import::WorkFlow::Helpers') or die;
my $helpers = Genome::InstrumentData::Command::Import::WorkFlow::Helpers->get;

my $tmp_dir = File::Temp::tempdir(CLEANUP => 1);

my @source_fastq_base_names = (qw/ input.1.fastq.gz input.2.fastq /);
my @source_fastq_paths = map { $tmp_dir.'/'.$_ } @source_fastq_base_names;
for my $source_fastq_base_name ( @source_fastq_base_names ) {
    my $source_fastq_path = $tmp_dir.'/'.$source_fastq_base_name;
    push @source_fastq_paths, $source_fastq_path;
    Genome::Sys->create_symlink($test_dir.'/'.$source_fastq_base_name, $source_fastq_path);
}

my $library = Genome::Library->__define__(
    id => -1, name => '__TEST_SAMPLE__-extlibs', sample => Genome::Sample->__define__(name => '__TEST_SAMPLE__')
);
ok($library, 'define library');

my $cmd = Genome::InstrumentData::Command::Import::WorkFlow::FastqsToBam->execute(
    working_directory => $tmp_dir,
    fastq_paths => \@source_fastq_paths,
    library => $library,
);
ok($cmd->result, 'execute');
my $output_bam_path = $cmd->output_bam_path;
is($output_bam_path, $tmp_dir.'/__TEST_SAMPLE__.bam', 'bam path named correctly');
ok(-s $output_bam_path, 'bam path exists');
is(File::Compare::compare($output_bam_path, $test_dir.'/input.fastq.unsorted.bam'), 0, 'bam matches');

for ( my $i = 0; $i < @source_fastq_paths; $i++ ) {
    ok(!glob($source_fastq_paths[$i].'*'), 'removed fastq '.($i + 1).' after conversion to bam');
}

#print "$tmp_dir\n"; <STDIN>;
done_testing();
