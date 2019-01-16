#!/usr/bin/env genome-perl

use strict;
use warnings;

BEGIN {
    $ENV{UR_DBI_NO_COMMIT} = 1;
};

use above 'Genome';
use Genome::Test::Factory::AnalysisProject;
use Test::More tests => 10;
use Genome::Utility::Test qw(compare_ok);

my $class = 'Genome::FeatureList::Command::DumpIntervalList';

use_ok($class);

Genome::Test::Factory::AnalysisProject->setup_system_analysis_project;
my $reference= Genome::Model::Build::ReferenceSequence->get(name => 'GRCh37-lite-build37');

my $bed_path = Genome::Sys->create_temp_file_path;
Genome::Sys->write_file($bed_path,
    "1\t11\t12\tregion1\n",
    "1\t40\t50\tregion2\n",
);
my $fl_cmd = Genome::FeatureList::Command::Create->create(
    name => 'test feature list for dump-interval-list',
    format => 'true-BED',
    file_path => $bed_path,
    content_type => 'targeted',
    description => 'just a test',
    source => 'GMS',
    reference => $reference,
);
my $test_fl = $fl_cmd->execute;
isa_ok($test_fl, 'Genome::FeatureList', 'setup test feature-list');

my $interval_list_path = Genome::Sys->create_temp_file_path;
my $dump_command = $class->create(
    feature_list => $test_fl,
    reference_build => $reference,
    track_name => 'target_region',
    output_path => $interval_list_path,
);
isa_ok($dump_command, $class, 'created interval-list command');
ok($dump_command->execute, 'executed interval-list command');

ok(-e $interval_list_path, 'created interval-list');

my $expected_file = Genome::Sys->create_temp_file_path;
Genome::Sys->write_file($expected_file,
"\@HD	VN:1.4	SO:unsorted\n",
"\@SQ etc. etc.\n" x 84,
"1	12	12	+	r0\n",
"1	41	50	+	r1\n"
);

compare_ok($interval_list_path, $expected_file, 
    name => 'produced expected interval-list',
    filters => [qr/\@SQ.*?$/],
);

my $interval_list_path2 = Genome::Sys->create_temp_file_path;
my $dump_command_with_preserved_regions = $class->create(
    feature_list => $test_fl,
    reference_build => $reference,
    track_name => 'target_region',
    output_path => $interval_list_path2,
    preserve_region_names => 1,
    merge => 0,
);
isa_ok($dump_command_with_preserved_regions, $class, 'created second interval-list command');
ok($dump_command_with_preserved_regions->execute, 'execute second interval-list command');

ok(-e $interval_list_path2, 'created interval-list');

my $expected_file2 = Genome::Sys->create_temp_file_path;
Genome::Sys->write_file($expected_file2,
"\@HD	VN:1.4	SO:unsorted\n",
"\@SQ etc. etc.\n" x 84,
"1	12	12	+	region1\n",
"1	41	50	+	region2\n"
);

compare_ok($interval_list_path2, $expected_file2,
    name => 'produced expected interval-list with region names',
    filters => [qr/\@SQ.*?$/],
);
