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

my $analysis_project = Genome::Config::AnalysisProject->__define__(name => 'TEST-AnP');
ok($analysis_project, 'define analysis project');
my $library = Genome::Library->__define__(name => 'TEST-sample-libs', sample => Genome::Sample->__define__(name => 'TEST-sample'));
ok($library, 'define library');
my @source_files = (qw/ in.1.fastq in.2.fastq /);
my %required_params = (
    analysis_project_id => $analysis_project->id,
    library_id => $library->id,
    source_paths => \@source_files,
);

my $inputs = $class->create(
    %required_params,
    instrument_data_properties => {
        description => 'imported',
        downsample_ratio => 0.7,
        import_source_name => 'TGI',
        this => 'that',
    },
);
ok($inputs, 'create inputs');
is($inputs->format, 'fastq', 'source files format is fastq');
is($inputs->library_name, $library->name, 'library_name');
is($inputs->sample_name, $library->sample->name, 'sample_name');

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

# instrument data
ok(!$inputs->instrument_data_for_original_data_path, 'no instrument_data_for_original_data_path ... yet');
my $instdata = Genome::InstrumentData::Imported->__define__;
ok($instdata, 'define instdata');
ok($instdata->original_data_path($inputs->source_files->original_data_path), 'add original_data_path');
is_deeply([$inputs->instrument_data_for_original_data_path], [$instdata], 'instrument_data_for_original_data_path');

# as_hashref
is_deeply(
    $inputs->as_hashref,
    {
        analysis_project => $analysis_project,
        downsample_ratio => $instrument_data_properties{downsample_ratio},
        instrument_data_properties => \%instrument_data_properties,
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
        sub { $class->create(%required_params); },
        qr/No $name given to work flow inputs\!/,
        "create failed w/o $name",
    );
    $required_params{$name} = $value;
}

done_testing();
