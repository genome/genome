#!/usr/bin/env genome-perl

use strict;
use warnings;

use above 'Genome';

use Genome::Utility::Test;
use File::Spec;
use Test::Exception;
use Test::More tests => 9;

my $class = 'Genome::InstrumentData::Command::Import::Inputs::Factory';
use_ok($class) or die;
use_ok('Genome::InstrumentData::Command::Import::Inputs') or die;

my $data_dir = Genome::Utility::Test->data_dir_ok('Genome::InstrumentData::Command::Import', 'generate-cmds');
my $input_file = File::Spec->join($data_dir, 'info.tsv');
my $analysis_project = Genome::Config::AnalysisProject->create(name => '__TEST_AP__');
my $process = Genome::InstrumentData::Command::Import::Process->create(
    analysis_project => $analysis_project,
    import_file => $input_file,
);

my $factory = $class->create(
    process => $process,
    file => $input_file,
);
ok($factory, 'create factory');

subtest 'next' => sub{
    plan tests => 4;

    $factory->_line_number(1); # set line number
    my $inputs = $factory->next; # line 2
    ok($inputs, 'inputs from next (2)');
    is($inputs->line_number, 2, 'line_number is correct');
    is_deeply([$inputs->source_paths], ['bam2.bam'], 'source_paths is correct');
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

my $inputs_from_line_4;
subtest 'from_line' => sub{
    plan tests => 4;

    $inputs_from_line_4 = $factory->from_line_number(4);
    ok($inputs_from_line_4, 'inputs from_line_number (4)');
    is($inputs_from_line_4->line_number, 4, 'line_number is correct');
    is_deeply([$inputs_from_line_4->source_paths], ['bam3.bam'], 'source_paths is correct');
    is_deeply(
        $inputs_from_line_4->entity_params,
        {
            individual => { name => 'TeSt-0000', nomenclature => 'TeSt', upn => '0000', },
            sample => { name => 'TeSt-0000-01', nomenclature => 'TeSt', },
            library => { name => 'TeSt-0000-01-extlibs', },
            instdata => { lane => 7, downsample_ratio => '.1', original_data_path => 'bam3.bam', process_id => $process->id, },
        },
        'entity_params are correct',
    );
};

subtest 'from_params' => sub{
    plan tests => 11;

    my $sample_name = 'TEST-01-001';
    my $library = Genome::Library->__define__(name => $sample_name.'-libs', sample => Genome::Sample->__define__(name => $sample_name));
    ok($library, 'define library');

    my @source_files = (qw/ in.1.fastq in.2.fastq /);
    my $instrument_data_properties = {
        description => 'imported',
        downsample_ratio => 0.7,
        import_source_name => 'TGI',
        this => 'that',
    };

    my $inputs = $factory->from_params({
            base_working_directory => '/tmp',
            source_paths => \@source_files,
            entity_params => {
                individual => { name => 'TEST-01', nomenclature => 'TEST', upn => '01', },
                sample => { name => $sample_name, nomenclature => 'TEST', },
                library => { name => $sample_name.'-libs', },
                instdata => $instrument_data_properties,
            },
        });
    ok($inputs, 'create inputs');
    is($inputs->format, 'fastq', 'source files format is fastq');
    is($inputs->library, $library, 'library');
    is($inputs->base_working_directory, '/tmp', 'base_working_directory');

    # instrument data
    ok(!$inputs->instrument_data_for_original_data_path, 'no instrument_data_for_original_data_path ... yet');
    my $instdata = Genome::InstrumentData::Imported->__define__;
    ok($instdata, 'define instdata');
    ok($instdata->original_data_path($inputs->source_files->original_data_path), 'add original_data_path');
    is_deeply([$inputs->instrument_data_for_original_data_path], [$instdata], 'instrument_data_for_original_data_path');

    # as_hashref
    $instrument_data_properties->{process_id} = $process->id;
    $instrument_data_properties->{original_data_path} = join(',', @source_files);
    is_deeply($inputs->instrument_data_properties, $instrument_data_properties, 'instrument_data_properties');
    is_deeply(
        $inputs->as_hashref,
        {
            analysis_project => $analysis_project,
            downsample_ratio => 0.7,
            instrument_data_properties => $instrument_data_properties,
            library => $library,
            source_paths => \@source_files,
        },
        'inputs as_hashref',
    );
};

subtest 'from_inputs_id' => sub {
    plan tests => 5;

    my $inputs_from_cache = Genome::InstrumentData::Command::Import::Inputs->get( join("\t", $process->id, 4) );
    is($inputs_from_cache, $inputs_from_line_4, 'from_inputs_id for process id and line 4');

    my $inputs_from_process_and_line_number = $factory->from_inputs_id( join("\t", $process->id, 3) );
    ok($inputs_from_process_and_line_number, 'inputs from_line_number (3)');
    is($inputs_from_process_and_line_number->line_number, 3, 'line_number is correct');
    is_deeply([$inputs_from_process_and_line_number->source_paths], ['bam3.bam'], 'source_paths is correct');
    is_deeply(
        $inputs_from_process_and_line_number->entity_params,
        {
            individual => { name => 'TeSt-0000', nomenclature => 'TeSt', upn => '0000', },
            sample => { name => 'TeSt-0000-01', nomenclature => 'TeSt', },
            library => { name => 'TeSt-0000-01-extlibs', },
            instdata => { lane => 7, downsample_ratio => '.25', original_data_path => 'bam3.bam', process_id => $process->id, },
        },
        'entity_params are correct',
    );
};

subtest 'fails' => sub{
    plan tests => 9;

    throws_ok(
        sub{ $class->create(process => $process, file => File::Spec->join($data_dir, 'samples.blah')); },
        qr/Cannot determine type for file: .+. It needs to end with \.csv or \.tsv\./,
        'failed w/ invalid file type',
    );

    throws_ok(
        sub{ $class->create(process => $process, file => File::Spec->join($data_dir, 'empty.csv')); },
        qr/File \(.+\) is empty\!/,
        'failed w/ empty file',
    );

    throws_ok(
        sub{ $class->create(process => $process, file => File::Spec->join($data_dir, 'invalid-entity-type.csv')); },
        qr/Invalid entity type: unknown/,
        'failed w/ invalid entity type',
    );

    throws_ok(
        sub{ $class->create(process => $process, file => File::Spec->join($data_dir, 'invalid-sample-name.csv'))->next; },
        qr/Invalid sample name: INVALID.NAME-. It must have at least 3 parts separated by dashes./,
        'failed w/ invalid sample name',
    );

    throws_ok(
        sub{ $class->create(process => $process, file => File::Spec->join($data_dir, 'no-nomenclature.csv'))->next; },
        qr/No sample\.nomenclature column given\! It is required to resolve entity names when no sample name is given\./,
        'failed w/o sample name and no sample.nomenclature',
    );

    throws_ok(
        sub{ $class->create(process => $process, file => File::Spec->join($data_dir, 'no-individual-name-part.csv'))->next; },
        qr/No individual\.name_part column_given! It is required to resolve entity names when no sample name is given\./,
        'failed w/o sample name and no individual.name_part',
    );

    throws_ok(
        sub{
            $class->create(process => $process, file => File::Spec->join($data_dir, 'no-sample-name-part.csv'))->next; },
        qr/No sample\.name_part column_given! It is required to resolve entity names when no sample name is given\./,
        'failed w/o sample name and no sample.name_part',
    );

    throws_ok(
        sub{ $class->create(process => $process, file => File::Spec->join($data_dir, 'individual-name-mismatch.csv'))->next; },
        qr/Invalid individual name: TGI-AAAA\. It must include the first part of the sample name: TGI-AA12345-Z98765\./,
        'failed when sample name does not include individual name',
    );

    throws_ok(
        sub{ $class->create->from_params({ base_working_directory => '/dev/null', entity_params => {},}); },
        qr/Can't validate/,
        'failed when given non writable base_working_directory',
    );

};

done_testing();
