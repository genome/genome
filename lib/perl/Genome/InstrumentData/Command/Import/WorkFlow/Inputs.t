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

my $inputs = $class->create(
    source_files => 'in.1.fastq,in.2.fastq',
    instrument_data_properties => [qw/ 
        description=imported
        downsample_ratio=0.7
        import_source_name=TGI
        this=that
    /],
);
ok($inputs, 'create inputs');
is_deeply(
    $inputs->for_worklflow,
    { 
        source_files => [qw/ in.1.fastq in.2.fastq /],
        instrument_data_properties => {
            downsample_ratio => 0.7,
            description => 'imported',
            import_source_name => 'TGI',
            original_data_path => 'in.1.fastq,in.2.fastq',
            this => 'that', 
        },
    },
    'inputs for workflow',
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
            source_files => 'in.bam',
            instrument_data_properties => [qw/ 
            description=imported
            description=inported
            /],
        );
    },
    qr/Multiple values for instrument data property! description => imported, inported/,
    "execute failed w/ duplicate key, diff value",
);

throws_ok(
    sub{
        $class->create(
            source_files => 'in.bam',
            instrument_data_properties => [qw/ 
            description=
            /],
        );
    },
    qr#Failed to parse with instrument data property label/value! description=#,
    "execute failes w/ missing value",
);

done_testing();
