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
use LWP::Simple;
use Sub::Install;
use Test::Exception;
use Test::MockObject;
use Test::More;

my $class = 'Genome::InstrumentData::Command::Import::WorkFlow::SourceFile';
use_ok($class) or die;
my $test_dir = Genome::Utility::Test->data_dir_ok('Genome::InstrumentData::Command::Import') or die;

for my $format (qw/ bam fastq sra /) {
    my $source_file = 'file.'.$format;
    my $sf = $class->create(path => $source_file);
    ok($sf, 'create');
    is($sf->format, $format, 'format');
    is($sf->retrieval_method, 'local disk', 'retrieval_method');
    is($sf->md5_path, $source_file.'.md5', 'md5_path');
    is($sf->original_md5_path, $source_file.'.md5-orig', 'original_data_path_md5');
}

# Remote bam
my $headers = Test::MockObject->new;
$headers->set_always('content_length', 2048);
my $response = Test::MockObject->new;
$response->set_true('is_success');
$response->set_always('headers', $headers);
my $lwp = Test::MockObject->new;
$lwp->set_always('head', $response);
Sub::Install::reinstall_sub({
        code => sub{ return $lwp; },
        into => 'LWP::UserAgent',
        as => 'new',
    });

my $sf = $class->create(path => 'http://file.bam');
ok($sf, 'create');
is($sf->format, 'bam', 'format for remote bam');
is($sf->retrieval_method, 'remote url', 'retrieval_method for remote bam');
is($sf->file_size, 2048, 'file_size');
is($sf->kilobytes_required_for_processing, 9, 'kilobytes_required_for_processing for bam');

# ERRORS
throws_ok(
    sub{ $class->create; },
    qr/No source file given\!/,
    'create fails w/o source file',
);

done_testing();
