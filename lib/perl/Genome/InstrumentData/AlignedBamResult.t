#! /gsc/bin/perl

use strict;
use warnings;

use above 'Genome';

require File::Temp;
require Genome::Utility::Test;
use Test::More;

my $class = 'Genome::InstrumentData::AlignedBamResult';
use_ok($class) or die;
my $data_dir = Genome::Utility::Test->data_dir_ok($class);
class Genome::InstrumentData::AlignedBamResultTester {
    is => $class,
};

my $tmpdir = File::Temp::tempdir(CLEANUP => 1);
my $aligned_bam_result = Genome::InstrumentData::AlignedBamResultTester->__define__(
    id => -1337,
    output_dir => $tmpdir,
);
ok($aligned_bam_result, 'defined aligned bam result');

my $bam_path = $aligned_bam_result->bam_path;
is($bam_path, $tmpdir.'/'.$aligned_bam_result->id.'.bam', 'bam path named correctly');
is($aligned_bam_result->bam_file, $bam_path, 'bam file is same as bam path');
Genome::Sys->create_symlink($data_dir.'/input.bam', $bam_path);
ok(-s $bam_path, 'linked bam');

done_testing();
