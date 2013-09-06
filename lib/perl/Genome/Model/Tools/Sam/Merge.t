#!/usr/bin/env genome-perl

use strict;
use warnings;

use above "Genome";
use Genome::Model::Tools::Sam::Merge;
use Test::More;
use File::Compare;
use File::Copy;

if (`uname -a` =~ /x86_64/){
    plan tests => 9;
} else{
    plan skip_all => 'Must run on a 64 bit machine';
}

my $dir = $ENV{GENOME_TEST_INPUTS} . '/Genome-Tools-Sam-Merge';

my $input_normal = $dir. '/normal.tiny.bam';
my $input_tumor  = $dir. '/tumor.tiny.bam';
my $bam_index    = $dir. '/normal_tumor.tiny.bam.bai';

# step 1: test 1 file case

my $out_1_dir = File::Temp::tempdir(
    "SamMerge1_XXXXXX",
    TMPDIR => 1,
    CLEANUP => 1,
);

my $merged_file1 = $out_1_dir.'/merged_1.bam';

my $cmd_1 = Genome::Model::Tools::Sam::Merge->create(
    files_to_merge => [$input_normal],
    merged_file    => $merged_file1,
    bam_index      => 0,
    merger_name    => 'picard',
    merger_version => '1.22',
    use_version    => 'r544'
);

ok($cmd_1, "created command",);
ok($cmd_1->execute, "executed");
ok(-s $merged_file1, "output file is nonzero");
ok(!-s $merged_file1.'.bai', 'Turn off .bai bam index generation');

# step 1: test >1 input file case

my $out_2_dir = File::Temp::tempdir(
    "SamMerge2_XXXXXX",
    TMPDIR => 1,
    CLEANUP => 1,
);
my $merged_file2 = $out_2_dir.'/merged_2.bam';

my $input_normal2 = $out_2_dir.'/normal.tiny.bam';
my $input_tumor2  = $out_2_dir.'/tumor.tiny.bam';
copy $input_normal, $input_normal2;
copy $input_tumor, $input_tumor2;

ok(-s $input_normal2 and -s $input_tumor2, 'normal, tumor bam are copied to tmp dir ok');

my $cmd_2 = Genome::Model::Tools::Sam::Merge->create(
    files_to_merge => [$input_normal2, $input_tumor2],
    merged_file    => $merged_file2,
    merger_name    => 'picard',
    merger_version => '1.22',
    use_version    => 'r544'
);

ok($cmd_2, "created command");
ok($cmd_2->execute, "executed");
ok(-s $merged_file2, "output file is nonzero");
is(compare($bam_index, $merged_file2.'.bai'), 0, 'The bam index is generated as expected'); 
