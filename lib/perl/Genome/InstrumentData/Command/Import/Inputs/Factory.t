#!/usr/bin/env genome-perl

use strict;
use warnings;

use above 'Genome';

use Genome::Utility::Test;
use File::Spec;
use Test::Exception;
use Test::More tests => 6;

my $class = 'Genome::InstrumentData::Command::Import::Inputs::Factory';
use_ok($class) or die;

my $data_dir = Genome::Utility::Test->data_dir_ok('Genome::InstrumentData::Command::Import', 'generate-cmds');
my $input_file = File::Spec->join($data_dir, 'info.tsv');

my $factory;
subtest 'create and set_file' => sub{
    plan tests => 2;

    $factory = $class->create;
    ok($factory, 'create inputs factory');
    ok($factory->set_file($input_file), 'set_file');

};

my $expected_line4_ref = {
    file => $input_file,
    line_number => 4,
    individual => { name => 'TeSt-0000', nomenclature => 'TeSt', upn => '0000', },
    sample => { name => 'TeSt-0000-01', nomenclature => 'TeSt', },
    library => { name => 'TeSt-0000-01-extlibs', },
    instdata => {
        lane => 7,
        downsample_ratio => '.1',
    },
    source_files => ['bam3.bam'],
};

subtest 'next' => sub{
    plan tests => 2;

    for (1..3) { $factory->next }
    my $ref = $factory->next;
    is_deeply(
        $ref,
        $expected_line4_ref,
        'last ref is correct',
    );
    ok(!$factory->next, 'reached end of file');
};

subtest 'from_line' => sub{
    plan tests => 1;

    is_deeply(
        $factory->from_line_number(4),
        $expected_line4_ref,
        'from_line 4 return is correct',
    );
};


subtest 'fails' => sub{
    plan tests => 8;

    throws_ok(
        sub{ $factory->set_file(File::Spec->join($data_dir, 'samples.blah')); },
        qr/Cannot determine type for file: .+. It needs to end with \.csv or \.tsv\./,
        'failed w/ invalid file type',
    );

    throws_ok(
        sub{ $factory->set_file(File::Spec->join($data_dir, 'empty.csv')); },
        qr/File \(.+\) is empty\!/,
        'failed w/ empty file',
    );

    throws_ok(
        sub{ $factory->set_file(File::Spec->join($data_dir, 'invalid-entity-type.csv')); },
        qr/Invalid entity type: unknown/,
        'failed w/ invalid entity type',
    );

    throws_ok(
        sub{ 
            my $f = $class->create;
            $f->set_file(File::Spec->join($data_dir, 'invalid-sample-name.csv'));
            $f->next;
        },
        qr/Invalid sample name: INVALID.NAME-. It must have at least 3 parts separated by dashes./,
        'failed w/ invalid sample name',
    );

    throws_ok(
        sub{
            my $f = $class->create;
            $f->set_file(File::Spec->join($data_dir, 'no-nomenclature.csv'));
            $f->next;
        },
        qr/No sample\.nomenclature column given\! It is required to resolve entity names when no sample name is given\./,
        'failed w/o sample name and no sample.nomenclature',
    );

    throws_ok(
        sub{
            my $f = $class->create;
            $f->set_file(File::Spec->join($data_dir, 'no-individual-name-part.csv'));
            $f->next;
        },
        qr/No individual\.name_part column_given! It is required to resolve entity names when no sample name is given\./,
        'failed w/o sample name and no individual.name_part',
    );

    throws_ok(
        sub{
            my $f = $class->create;
            $f->set_file(File::Spec->join($data_dir, 'no-sample-name-part.csv'));
            $f->next;
        },
        qr/No sample\.name_part column_given! It is required to resolve entity names when no sample name is given\./,
        'failed w/o sample name and no sample.name_part',
    );

    throws_ok(
        sub{ 
            my $f = $class->create;
            $f->set_file(File::Spec->join($data_dir, 'individual-name-mismatch.csv'));
            $f->next;
        },
        qr/Invalid individual name: TGI-AAAA\. It must include the first part of the sample name: TGI-AA12345-Z98765\./,
        'failed when sample name does not include individual name',
    );
};

done_testing();
