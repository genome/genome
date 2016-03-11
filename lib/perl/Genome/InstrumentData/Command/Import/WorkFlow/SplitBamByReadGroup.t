#!/usr/bin/env genome-perl

use strict;
use warnings;

use above 'Genome';

require Genome::Utility::Test;
require File::Compare;
require File::Spec;
require File::Temp;
require List::MoreUtils;
require Sub::Install;
use Test::More tests => 25;

my $class = 'Genome::InstrumentData::Command::Import::WorkFlow::SplitBamByReadGroup';
use_ok($class) or die;
my $test_dir = Genome::Utility::Test->data_dir_ok('Genome::InstrumentData::Command::Import', 'v2') or die;
Genome::InstrumentData::Command::Import::WorkFlow::Helpers->overload_uuid_generator_for_class($class);

subtest 'read1 or read2' => sub{
    plan tests => 3;

    is(Genome::InstrumentData::Command::Import::WorkFlow::SplitBamByReadGroup::_read1_or_read2(77), 'read1', 'flag marked as read 1 is read1');
    is(Genome::InstrumentData::Command::Import::WorkFlow::SplitBamByReadGroup::_read1_or_read2(0), 'read1', 'unmarked flag is read1');
    is(Genome::InstrumentData::Command::Import::WorkFlow::SplitBamByReadGroup::_read1_or_read2(141), 'read2', 'flag marked as read 2 is read2');

};

my $tmp_dir = File::Temp::tempdir(CLEANUP => 1);
my $multi_rg_base_name = 'input.rg-multi.bam';
my $multi_rg_bam_path = File::Spec->join($tmp_dir, $multi_rg_base_name);
Genome::Sys->create_symlink(File::Spec->join($test_dir, $multi_rg_base_name), $multi_rg_bam_path);
ok(-s $multi_rg_bam_path, 'linked two read groups bam');
my $cmd = Genome::InstrumentData::Command::Import::WorkFlow::SplitBamByReadGroup->execute(
    working_directory => $tmp_dir,
    bam_path => $multi_rg_bam_path,
);
ok($cmd->result, 'execute');

my @output_bam_paths = $cmd->output_bam_paths;
my @bam_basenames = (qw/ 2883581797.paired 2883581797.read1 2883581797.read2 2883581798.paired 2883581798.read1 2883581798.read2 /);
is(@output_bam_paths, @bam_basenames, '6 read group bam paths');
for my $basename ( @bam_basenames ) {
    my $output_bam_path = File::Spec::->join($tmp_dir, 'input.rg-multi.'.$basename.'.bam');
    ok((List::MoreUtils::any { $_ eq $output_bam_path } @output_bam_paths), 'expected bam in output bams');
    ok(-s $output_bam_path, 'expected bam path exists');
    my $expected_bam_path = File::Spec->join($test_dir, 'split-by-rg.'.$basename.'.bam');
    is(File::Compare::compare($output_bam_path, $expected_bam_path), 0, "expected $basename bam path matches");
}
ok(!glob($multi_rg_bam_path.'*'), 'removed bam path and auxiliary files after spliting');

#print "$tmp_dir\n"; <STDIN>;
done_testing();
