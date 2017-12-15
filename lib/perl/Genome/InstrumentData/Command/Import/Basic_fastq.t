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
Genome::InstrumentData::Command::Import::WorkFlow::Helpers->overload_uuid_generator_for_class('Genome::InstrumentData::Command::Import::WorkFlow::FastqsToBam');

my $test_dir = Genome::Utility::Test->data_dir_ok('Genome::InstrumentData::Command::Import', 'v5');
my @source_files = (
    $test_dir.'/input.1.fastq.gz', 
    $test_dir.'/input.2.fastq',
);
is(grep({ -s $_ } @source_files), 2, 'source fastqs exist');

my $analysis_project = Genome::Config::AnalysisProject->create(name => '__TEST_AP__');
ok($analysis_project, 'create analysis project');
my $library = Genome::Library->create(
    name => '__TEST_SAMPLE__-extlibs', sample => Genome::Sample->create(name => '__TEST_SAMPLE__')
);
ok($library, 'Create library');

my $wd = Genome::Sys->create_temp_directory();

my $cmd = Genome::InstrumentData::Command::Import::Basic->create(
    analysis_project => $analysis_project,
    library => $library,
    source_files => \@source_files,
    import_source_name => 'broad',
    instrument_data_properties => [qw/ lane=2 flow_cell_id=XXXXXX /],
    base_working_directory => $wd,
);
ok($cmd, "create import command");
ok($cmd->execute, "execute import command");

my %instrument_data;
for my $source_file ( @source_files ) {
    my $md5 = Genome::InstrumentData::Command::Import::WorkFlow::Helpers->load_md5($source_file.'.md5');
    ok($md5, 'load source md5');
    my @inst_data = map { $_->instrument_data } Genome::InstrumentDataAttribute->get(
        attribute_label => 'original_data_path_md5',
        attribute_value => $md5,
    );
    is(@inst_data, 1, "got instrument data for md5 $md5") or die;
    for ( @inst_data ) { $instrument_data{$_->id} = $_; }
}
my @instrument_data = values %instrument_data;
is(@instrument_data, 1, "got instrument data for original fastq md5s") or die;
my $instrument_data = $instrument_data[0];
is($instrument_data->original_data_path, join(',', @source_files), 'original_data_path correctly set');
is($instrument_data->import_format, 'bam', 'import_format is bam');
is($instrument_data->sequencing_platform, 'solexa', 'sequencing_platform correctly set');
is($instrument_data->is_paired_end, 1, 'is_paired_end correctly set');
is($instrument_data->read_count, 2000, 'read_count correctly set');
is($instrument_data->read_length, 75, 'read_length correctly set');
is(eval{ $instrument_data->attributes(attribute_label => 'lane')->attribute_value }, 2, 'lane correctly set');
is(eval{ $instrument_data->attributes(attribute_label => 'flow_cell_id')->attribute_value }, 'XXXXXX', 'flow_cell_id correctly set');
is($instrument_data->analysis_projects, $analysis_project, 'set analysis project');

my $bam_path = $instrument_data->bam_path;
ok(-s $bam_path, 'bam path exists');
is($bam_path, $instrument_data->data_directory.'/all_sequences.bam', 'bam path correctly named');
is(eval{$instrument_data->attributes(attribute_label => 'bam_path')->attribute_value}, $bam_path, 'set attributes bam path');

# Rely on flagstat until a better BAM comparison process is defined
# is(File::Compare::compare($bam_path, $test_dir.'/all_sequences.basic-fastq-t.bam'), 0, 'bam matches');

is(File::Compare::compare($bam_path.'.flagstat', $test_dir.'/all_sequences.basic-fastq-t.bam.flagstat'), 0, 'flagstat matches');

my $allocation = $instrument_data->disk_allocation;
ok($allocation, 'got allocation');
ok($allocation->kilobytes_requested > 0, 'allocation kb was set');

#print $instrument_data->data_directory."\n";<STDIN>;
done_testing();
