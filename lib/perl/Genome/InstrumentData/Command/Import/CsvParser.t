#!/usr/bin/env genome-perl

use strict;
use warnings;

use above 'Genome';

use Genome::Utility::Test;
use File::Spec;
use Test::Exception;
use Test::More;

my $class = 'Genome::InstrumentData::Command::Import::CsvParser';
use_ok($class) or die;

my $data_dir = Genome::Utility::Test->data_dir_ok('Genome::InstrumentData::Command::Import', 'generate-cmds');
my $input_file = File::Spec->join($data_dir, 'info.csv');

my $csv_parser = $class->create(file => File::Spec->join($data_dir, 'info.tsv'));
my $cnt = 0;
my $ref;
while ( $ref = $csv_parser->next ) {
    $cnt++;
}
is($cnt, 4, 'read 4 entries');
is_deeply(
    $ref,
    {
        'library.name' => 'TeSt-0000-01-extlibs',
        'instdata.source_files' => 'bam3.bam',
        'instdata.lane' => 7,
        'instdata.downsample_ratio' => '.1',
    }, 
    'last ref is correct',
);

# Fails
## invalid file type
throws_ok(
    sub{ $class->create(file => File::Spec->join($data_dir, 'samples.blah')); },
    qr/Cannot determine type for file: .+. It needs to end with \.csv or \.tsv\./,
    'failed w/ invalid file type',
);

## empty file
throws_ok(
    sub{ $class->create(file => File::Spec->join($data_dir, 'empty.csv')); },
    qr/File \(.+\) is empty\!/,
    'failed w/ empty file',
);

## invalid entity type
throws_ok(
    sub{ $class->create(file => File::Spec->join($data_dir, 'invalid-entity-type.csv')); },
    qr/Invalid entity type: unknown/,
    'failed w/ invalid entity type',
);

## invalid sample name
throws_ok(
    sub{ $class->create(file => File::Spec->join($data_dir, 'invalid-sample-name.csv'))->next; },
    qr/Invalid sample name: INVALID.NAME-. It must have at least 3 parts separated by dashes./,
    'failed w/ invalid sample name',
);

## no sample name then nomenclature is required
throws_ok(
    sub{ $class->create(file => File::Spec->join($data_dir, 'no-nomenclature.csv')); },
    qr/No sample\.nomenclature column given\! It is required to resolve entity names when no sample name is given\./,
    'failed w/o sample name and sample.nomenclature',
);

done_testing(); exit;

## no sample name then individual name part is required
throws_ok(
    sub{ $class->create(file => File::Spec->join($data_dir, 'no-individual-name-part.csv')); },
    qr/No individual\.name_part column_given! It is required to resolve entity names when no sample name is given\./,
    'failed w/o sample name and individual.name_part',
);

## no sample name then sample name part is required
throws_ok(
    sub{ $class->create(file => File::Spec->join($data_dir, 'no-sample-name-part.csv')); },
    qr/No sample\.name_part column_given! It is required to resolve entity names when no sample name is given\./,
    'failed w/o sample name and sample.name_part',
);

## no sample name then sample name part is required
throws_ok(
    sub{ $class->create(file => File::Spec->join($data_dir, 'individual-name-mismatch.csv')); },
    qr/Invalid individual name: TGI-AAAA\. It must include the first part of the sample name: TGI-AA12345-Z98765\./,
    'failed when sample name does not include individual name',
);

## missing required property
throws_ok(
    sub{ $class->create(file => File::Spec->join($data_dir, 'missing-required-property.csv')); },
    qr/Missing required individual properties\: taxon/,
    'failed when missing requried property [individual taxon]',
);

## unknown property
throws_ok(
    sub{ $class->create(file => File::Spec->join($data_dir, 'unknown-property.csv')); },
    qr/Unknown individual property\: unknown_property/,
    'failed w/ unknown property',
);

done_testing();
