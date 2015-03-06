#!/usr/bin/env genome-perl

BEGIN {
    $ENV{UR_DBI_NO_COMMIT} = 1;
}

use strict;
use warnings;

use above 'Genome';
require File::Spec;
use Genome::Utility::Test;
use Test::Exception;
use Test::More;

use_ok('Genome::InstrumentData::Command::Import::GenerateSourceFilesTsv') or die;
my $data_dir = Genome::Utility::Test->data_dir_ok('Genome::InstrumentData::Command::Import', 'generate-cmds');
my $input_file = File::Spec->join($data_dir, 'info.csv');

my $tmpdir = Genome::Sys->create_temp_directory;
my $source_files_tsv_name = 'source-files.tsv';
my $source_files_tsv = File::Spec->join($tmpdir, $source_files_tsv_name);
my $cmd = Genome::InstrumentData::Command::Import::GenerateSourceFilesTsv->execute(
    file => $input_file,
    output_file => $source_files_tsv,
);
ok($cmd->result, 'execute');
Genome::Utility::Test::compare_ok($cmd->output_file, File::Spec->join($data_dir, $source_files_tsv_name), 'source-files.tsv files matche');

# Fails 
## no source files in header
throws_ok(
    sub{ Genome::InstrumentData::Command::Import::GenerateSourceFilesTsv->execute(file => File::Spec->join($data_dir, 'no-source-files-header.csv')); },
    qr/No source_files attribute for instdata in file/,
    'failed w/o source files in header',
);

## no source files
throws_ok(
    sub{ Genome::InstrumentData::Command::Import::GenerateSourceFilesTsv->execute(file => File::Spec->join($data_dir, 'no-source-files.csv')); },
    qr/No instrument data source files speicifed\!/,
    'failed w/o source files',
);

done_testing();
