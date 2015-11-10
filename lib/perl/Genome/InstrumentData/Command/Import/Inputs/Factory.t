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
my $analysis_project = Genome::Config::AnalysisProject->create(name => '__TEST_AP__');
my $process = Genome::InstrumentData::Command::Import::Process->create(
    analysis_project => $analysis_project,
    import_file => $input_file,
);

my $factory;
subtest 'create and set_file' => sub{
    plan tests => 2;

    $factory = $class->create(process => $process);
    ok($factory, 'create inputs factory');
    ok($factory->set_file($input_file), 'set_file');

};

subtest 'next' => sub{
    plan tests => 4;

    $factory->_line_number(1); # set line number
    my $inputs = $factory->next; # line 2
    ok($inputs, 'inputs from next (2)');
    is($inputs->line_number, 2, 'line_number is correct');
    is_deeply($inputs->source_paths, ['bam2.bam'], 'source_paths is correct');
    is_deeply(
        $inputs->entity_params,
        {
            individual => { name => 'TeSt-0000', nomenclature => 'TeSt', upn => '0000', },
            sample     => { name => 'TeSt-0000-00', nomenclature => 'TeSt', },
            library    => { name => 'TeSt-0000-00-extlibs', },
            instdata   => { lane => 8, process_id => $process->id, original_data_path => 'bam2.bam', },
        },
        'entity_params are correct',
    );

};

subtest 'from_line' => sub{
    plan tests => 4;

    my $inputs = $factory->from_line_number(4);
    ok($inputs, 'inputs from_line_number (4)');
    is($inputs->line_number, 4, 'line_number is correct');
    is_deeply($inputs->source_paths, ['bam3.bam'], 'source_paths is correct');
    is_deeply(
        $inputs->entity_params,
        {
            individual => { name => 'TeSt-0000', nomenclature => 'TeSt', upn => '0000', },
            sample => { name => 'TeSt-0000-01', nomenclature => 'TeSt', },
            library => { name => 'TeSt-0000-01-extlibs', },
            instdata => { lane => 7, downsample_ratio => '.1', original_data_path => 'bam3.bam', process_id => $process->id, },
        },
        'entity_params are correct',
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
