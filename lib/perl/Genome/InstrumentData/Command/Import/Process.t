#!/usr/bin/env genome-perl

BEGIN {
    $ENV{UR_DBI_NO_COMMIT} = 1;
    $ENV{UR_COMMAND_DUMP_STATUS_MESSAGES} = 1;
    $ENV{UR_COMMAND_DUMP_DEBUG_MESSAGES} = 1;
}

use strict;
use warnings;

use above "Genome";

require Digest::MD5;
require File::Spec;
require Genome::Utility::Test;
use Test::More;
use Test::Exception;

my $class = 'Genome::InstrumentData::Command::Import::Process';
use_ok($class) or die;

my $data_dir = Genome::Utility::Test->data_dir_ok('Genome::InstrumentData::Command::Import', 'generate-cmds');
my $import_file = File::Spec->join($data_dir, 'info.tsv');
my $anp = Genome::Config::AnalysisProject->__define__(name => '__TEST_ANP__');

throws_ok(sub{ $class->create(analysis_project => $anp); }, qr/^No import file given to create process\!/, 'create w/o file');
throws_ok(sub{ $class->create(analysis_project => $anp, import_file => 'blah'); }, qr#Import file \(.+/blah\) given to create process does not exist\!#, 'create w/ file that does not exist');

my $process = $class->create(
    analysis_project => $anp,
    import_file => $import_file,
);
ok($process, 'create process');
is($process->analysis_project, $anp, 'analysis_project');
is($process->import_file, $import_file, 'import_file');
ok($process->import_md5, 'import_md5');
ok(!$process->saved_import_file, 'saved_import_file not defined b/c no allocation');
ok($process->create_disk_allocation, 'create_disk_allocation');
ok(-s $process->saved_import_file, 'saved_import_file');
my $instdata = Genome::InstrumentData::Imported->create(process_id => $process->id);
is_deeply([$process->instrument_data], [$instdata], 'instrument_data');

done_testing();
