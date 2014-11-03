#! /gsc/bin/perl

use strict;
use warnings;

use above 'Genome';

require Genome::Utility::Test;
require File::Compare;
require File::Spec;
require File::Temp;
use Test::More;

use_ok('Genome::InstrumentData::Command::Import::WorkFlow::SanitizeBam') or die;
my $test_dir = Genome::Utility::Test->data_dir_ok('Genome::InstrumentData::Command::Import', 'bam/v4') or die;

my $tmp_dir = File::Temp::tempdir(CLEANUP => 1);
my $dirty_bam_base_name = 'test.bam';
my $dirty_bam_path = File::Spec->catfile($tmp_dir, $dirty_bam_base_name);
Genome::Sys->create_symlink( File::Spec->catfile($test_dir, $dirty_bam_base_name), $dirty_bam_path );
ok(-s $dirty_bam_path, 'linked dirty bam path');

my $dirty_flagstat_base_name = $dirty_bam_base_name.'.flagstat';
my $dirty_flagstat_path = File::Spec->catfile($tmp_dir, $dirty_flagstat_base_name);
Genome::Sys->create_symlink( File::Spec->catfile($test_dir, $dirty_flagstat_base_name), $dirty_flagstat_path);
ok(-s $dirty_flagstat_path, 'linked dirty bam flagstat path');

my $cmd = Genome::InstrumentData::Command::Import::WorkFlow::SanitizeBam->execute(
    bam_path => $dirty_bam_path,
);
ok($cmd->result, 'execute');

my $clean_bam_base_name = 'test.clean.bam';
my $clean_bam_path = File::Spec->catfile($tmp_dir, $clean_bam_base_name);
is($cmd->output_bam_path, $clean_bam_path, 'clean bam path named correctly');
ok(-s $clean_bam_path, 'clean bam path exists');
my $expected_bam_path = File::Spec->catfile($test_dir, $clean_bam_base_name);
is(File::Compare::compare($clean_bam_path, $expected_bam_path), 0, 'clean bam matches');

my $clean_flagstat_base_name = $clean_bam_base_name.'.flagstat';
my $clean_flagstat_path = File::Spec->catfile($tmp_dir, $clean_flagstat_base_name);
ok(-s $clean_flagstat_path, 'flagstat path exists');
my $expected_flagstat_path = File::Spec->catfile($test_dir, $clean_flagstat_base_name);
is(File::Compare::compare($clean_flagstat_path, $expected_flagstat_path), 0, 'flagstat matches');
ok(!glob($dirty_bam_path.'*'), 'removed dirty bam path and auxillary files after sanitizing');

#print "$tmp_dir\n"; <STDIN>;
done_testing();
