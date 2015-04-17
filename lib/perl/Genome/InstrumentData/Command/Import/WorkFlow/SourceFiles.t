#!/usr/bin/env genome-perl

BEGIN {
    $ENV{UR_COMMAND_DUMP_STATUS_MESSAGES} = 1;
}

use strict;
use warnings;

use above 'Genome';

require File::Spec;
require File::Temp;
require Genome::Utility::Test;
use Test::Exception;
use Test::More;

my $class = 'Genome::InstrumentData::Command::Import::WorkFlow::SourceFiles';
use_ok($class) or die;
my $test_dir = Genome::Utility::Test->data_dir_ok('Genome::InstrumentData::Command::Import') or die;

my @source_files = map { File::Spec->join($test_dir, 'fastq', 'v3', $_) } (qw/ input.1.fastq.gz input.2.fastq /);

my $sf = $class->create(source_files => \@source_files);
ok($sf, 'create');
is($sf->format, 'fastq', 'format');
is($sf->retrieval_method, 'local disk', 'retrieval_method');
is($sf->kilobytes_required_for_processing, 746, 'kilobytes_required_for_processing');
ok($sf->verify_adequate_disk_space_is_available_for_processing('/tmp'), 'verify_adequate_disk_space_is_available_for_processing');

# ERRORS
throws_ok(
    sub{ $class->create(); },
    qr/No source files given\!/,
    'create fails w/o source files',
);

throws_ok(
    sub{ $class->create(source_files => [qw/ file.bam file.fastq /]); },
    qr/Mixed values for source file format\!/,
    'create fails w/ different formats',
);

throws_ok(
    sub{ $class->create(source_files => [qw{ http://file.bam file.bam }]); },
    qr/Mixed values for source file retrieval_method\!/,
    'create fails w/ different retrieval_methods',
);

done_testing();
