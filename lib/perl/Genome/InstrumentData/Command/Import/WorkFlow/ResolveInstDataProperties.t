#!/usr/bin/env genome-perl

use strict;
use warnings;

BEGIN {
    $ENV{UR_DBI_NO_COMMIT} = 1;
    $ENV{UR_COMMAND_DUMP_STATUS_MESSAGES} = 1;
};

use above "Genome";

use Test::More;

use_ok('Genome::InstrumentData::Command::Import::WorkFlow::ResolveInstDataProperties') or die;

my $source = join(',', (qw/ file1 file2 /));
my $cmd = Genome::InstrumentData::Command::Import::WorkFlow::ResolveInstDataProperties->execute(
    source => $source,
    instrument_data_properties => [qw/ 
        description=imported
        downsample_ratio=0.7
        import_source_name=TGI
        this=that
    /],
);
ok($cmd->result,"execute succeeded");
is_deeply(
    $cmd->resolved_instrument_data_properties,
    { 
        downsample_ratio => 0.7,
        description => 'imported',
        import_source_name => 'TGI',
        original_data_path => $source,
        this => 'that', 
    },
    'resolved_instrument_data_properties',
);

# ERRORS
$cmd = Genome::InstrumentData::Command::Import::WorkFlow::ResolveInstDataProperties->execute(
    source => $source,
    instrument_data_properties => [qw/ 
        description=imported
        description=inported
    /],
);
ok(!$cmd->result,"execute failed w/ duplicate key, diff value");
is($cmd->error_message, 'Multiple values for instrument data property! description => imported, inported', 'correct error');

$cmd = Genome::InstrumentData::Command::Import::WorkFlow::ResolveInstDataProperties->execute(
    source => $source,
    instrument_data_properties => [qw/ 
        description=
    /],
);
ok(!$cmd->result,"execute failes w/ missing value");
is($cmd->error_message, 'Failed to parse with instrument data property label/value! description=', 'correct error');

done_testing();
