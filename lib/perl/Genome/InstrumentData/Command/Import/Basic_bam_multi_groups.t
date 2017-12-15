#!/usr/bin/env genome-perl

BEGIN {
    $ENV{UR_USE_DUMMY_AUTOGENERATED_IDS} = 1;
    $ENV{UR_DBI_NO_COMMIT} = 1;
    $ENV{UR_COMMAND_DUMP_STATUS_MESSAGES} = 1;
}

use strict;
use warnings;

use above "Genome";

require Genome::Utility::Test;
require File::Compare;
use Test::More;

use_ok('Genome::InstrumentData::Command::Import::Basic') or die;
use_ok('Genome::InstrumentData::Command::Import::WorkFlow::Helpers') or die;
Genome::InstrumentData::Command::Import::WorkFlow::Helpers->overload_uuid_generator_for_class('Genome::InstrumentData::Command::Import::WorkFlow::SanitizeAndSplitBam');

my $analysis_project = Genome::Config::AnalysisProject->create(name => '__TEST_AP__');
ok($analysis_project, 'create analysis project');
my $library = Genome::Library->create(
    name => '__TEST_SAMPLE__-extlibs', sample => Genome::Sample->create(name => '__TEST_SAMPLE__')
);
ok($library, 'Create library');

my $test_dir = Genome::Utility::Test->data_dir_ok('Genome::InstrumentData::Command::Import', 'v02');
my $source_bam = $test_dir.'/input.multi-rg.bam';
ok(-s $source_bam, 'source bam exists') or die;

my $wd = Genome::Sys->create_temp_directory();

my $cmd = Genome::InstrumentData::Command::Import::Basic->create(
    analysis_project => $analysis_project,
    library => $library,
    source_files => [$source_bam],
    import_source_name => 'broad',
    instrument_data_properties => [qw/ lane=2 flow_cell_id=XXXXXX /],
    base_working_directory => $wd,
);
ok($cmd, "create import command");
ok($cmd->execute, "execute import command");

my $md5 = Genome::InstrumentData::Command::Import::WorkFlow::Helpers->load_md5($source_bam.'.md5');
ok($md5, 'load source md5');
my @instdata_md5_attr = Genome::InstrumentDataAttribute->get(
    attribute_label => 'original_data_path_md5',
    attribute_value => $md5,
);
is(@instdata_md5_attr, 4, "got instrument data for md5 $md5");

my %instrument_data = map { $_->attribute_value => $_->instrument_data } Genome::InstrumentDataAttribute->get(
    attribute_label => 'segment_id',
    instrument_data_id => [ map { $_->instrument_data_id } @instdata_md5_attr ],
);
my %expected_read_groups = (
    # 2883581797
    '11111111111111111111111111111111' => [qw/ paired    128 /],
    '22222222222222222222222222222222' => [qw/ singleton   1 /],
    # 2883581798
    '33333333333333333333333333333333' => [qw/ paired    124 /],
    '44444444444444444444444444444444' => [qw/ singleton   3 /],
);
is_deeply([sort keys %instrument_data], [sort keys %expected_read_groups], 'got instrument data for md5 and read groups');

for my $rg_id ( sort keys %expected_read_groups ) {
    my $instrument_data = $instrument_data{$rg_id};
    ok($instrument_data, "got instrument data for rg_id: $rg_id");

    my $type = $expected_read_groups{$rg_id}->[0];
    my $is_paired = $type eq 'paired' ? 1 : 0;
    is($instrument_data->is_paired_end, $is_paired, "$rg_id correct is_paired_end");

    is($instrument_data->original_data_path, $source_bam, 'original_data_path correctly set');
    is($instrument_data->import_format, 'bam', 'import_format is bam');
    is($instrument_data->sequencing_platform, 'solexa', 'sequencing_platform correctly set');
    is($instrument_data->read_count, $expected_read_groups{$rg_id}->[1], 'read_count correctly set');
    is($instrument_data->read_length, 100, 'read_length correctly set');
    is($instrument_data->analysis_projects, $analysis_project, 'set analysis project');

    my $bam_path = $instrument_data->bam_path;
    ok(-s $bam_path, 'bam path exists');
    is($bam_path, $instrument_data->data_directory.'/all_sequences.bam', 'bam path correctly named');
    is(eval{$instrument_data->attributes(attribute_label => 'bam_path')->attribute_value}, $bam_path, 'set attributes bam path');
    my $bam_basename = join('.', 'basic-bam-multi-rg', $rg_id, $type, 'bam');
    is(File::Compare::compare($bam_path, $test_dir.'/'.$bam_basename), 0, 'bam matches');
    is(File::Compare::compare($bam_path.'.flagstat', $test_dir.'/'.$bam_basename.'.flagstat'), 0, 'flagstat matches');

    my $allocation = $instrument_data->disk_allocation;
    ok($allocation, 'got allocation');
    ok($allocation->kilobytes_requested > 0, 'allocation kb was set');

}

done_testing();
