#!/usr/bin/env genome-perl

use strict;
use warnings;

BEGIN {
    $ENV{UR_DBI_NO_COMMIT} = 1;
    $ENV{UR_COMMAND_DUMP_STATUS_MESSAGES} = 1;
};

use above "Genome";

use Test::Exception;
use Test::More;

my $class = 'Genome::InstrumentData::Command::Import::WorkFlow::Inputs';
use_ok($class) or die;

my @source_files = (qw/ in.1.fastq in.2.fastq /);
my $inputs = $class->create(
    source_files => \@source_files,
    instrument_data_properties => [qw/ 
        description=imported
        downsample_ratio=0.7
        import_source_name=TGI
        this=that
    /],
);
ok($inputs, 'create inputs');

isa_ok($inputs->source_files, 'Genome::InstrumentData::Command::Import::WorkFlow::SourceFiles', 'set _source_files');
is($inputs->format, 'fastq', 'source files format is fastq');

my %instrument_data_properties = (
    downsample_ratio => 0.7,
    description => 'imported',
    import_source_name => 'TGI',
    original_data_path => join(',', @source_files),
    this => 'that', 
);
is_deeply(
    $inputs->instrument_data_properties,
    \%instrument_data_properties,
    'instrument_data_properties',
);

# ERRORS
throws_ok(
    sub { $class->create(); },
    qr/No source files\!/,
    "create failed w/o source files",
);

throws_ok(
    sub {
        $class->create(
            source_files => [qw/ in.bam /],
            instrument_data_properties => [qw/ foo=bar foo=baz /],
        );
    },
    qr/Multiple values for instrument data property! foo => bar, baz/,
    "execute failed w/ duplicate key, diff value for instdata properties",
);

throws_ok(
    sub{
        $class->create(
            source_files => [qw/ in.bam /],
            instrument_data_properties => [qw/ 
            description=
            /],
        );
    },
    qr#Failed to parse with instrument data property label/value! description=#,
    "execute failes w/ missing value",
);

done_testing();
