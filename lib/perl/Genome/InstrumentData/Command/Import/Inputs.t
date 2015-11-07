#!/usr/bin/env genome-perl

use strict;
use warnings;

BEGIN {
    $ENV{UR_DBI_NO_COMMIT} = 1;
    $ENV{UR_COMMAND_DUMP_STATUS_MESSAGES} = 1;
};

use above "Genome";

require Genome::Utility::Test;
use Test::Exception;
use Test::More;

my $class = 'Genome::InstrumentData::Command::Import::Inputs';
use_ok($class) or die;

my $analysis_project = Genome::Config::AnalysisProject->__define__(name => 'TEST-AnP');
ok($analysis_project, 'define analysis project');
my $sample_name = 'TEST-01-001';
my $library = Genome::Library->__define__(name => $sample_name.'-libs', sample => Genome::Sample->__define__(name => $sample_name));
ok($library, 'define library');
my @source_files = (qw/ in.1.fastq in.2.fastq /);
my %required_params = (
    source_paths => \@source_files,
);

my $data_dir = Genome::Utility::Test->data_dir_ok('Genome::InstrumentData::Command::Import', 'generate-cmds');
my $import_file = File::Spec->join($data_dir, 'info.tsv');
my $process = Genome::InstrumentData::Command::Import::Process->create(
    analysis_project => $analysis_project,
    import_file => $import_file,
);
my $line_number = 0;
my $instrument_data_properties = {
    description => 'imported',
    downsample_ratio => 0.7,
    import_source_name => 'TGI',
    this => 'that',
};

my $inputs = $class->create(
    %required_params,
    process_id => $process->id,
    line_number => ++$line_number,
    entity_params => {
        individual => { name => 'TEST-01', nomenclature => 'TEST', upn => '01', },
        sample => { name => $sample_name, nomenclature => 'TEST', },
        library => { name => $sample_name.'-libs', },
        instdata => $instrument_data_properties,
    },
);
ok($inputs, 'create inputs');
is($inputs->format, 'fastq', 'source files format is fastq');
is($inputs->library, $library, 'library');
is($inputs->library_name, $library->name, 'library_name');
is($inputs->sample_name, $library->sample->name, 'sample_name');
is_deeply(
    $class->get( join("\t", $process->id, $line_number) ),
    $inputs,
    'get inputs',
);

# md5 props
is($inputs->lib_and_source_file_md5sum, '7f7d81', 'md5_properties');

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
        library_name => $library->name,
        sample_name => $library->sample->name,
        source_paths => \@source_files,
    },
    'inputs as_hashref',
);

# ERRORS
for my $name ( sort keys %required_params ) {
    my $value = delete $required_params{$name};
    throws_ok(
        sub { $class->create(
                process_id => $process->id,
                line_number => ++$line_number,
                %required_params,
            ); },
        qr/No $name given to work flow inputs\!/,
        "create failed w/o $name",
    );
    $required_params{$name} = $value;
}

done_testing();
