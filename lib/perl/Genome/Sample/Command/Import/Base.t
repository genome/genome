#!/usr/bin/env genome-perl

BEGIN {
    $ENV{UR_DBI_NO_COMMIT} = 1;
    $ENV{UR_COMMAND_DUMP_STATUS_MESSAGES} = 1;
}

use strict;
use warnings;

use above 'Genome';

use Test::More;

use_ok('Genome::Sample::Command::Import::Base') or die;
class SampleImporter {
    is => 'Genome::Sample::Command::Import::Base',
    has => {
        nomenclature => {},
    },
};

my $sample_name = '_SAMPLE_';
my $sample = Genome::Sample->__define__(
    name => $sample_name,
);
ok($sample, 'defined sample');

my $cmd = SampleImporter->execute(name => $sample_name, nomenclature => 'TGI');
ok($cmd->result, 'execute sample importer for existoing sample');
is($cmd->_sample, $sample, 'set sample');

done_testing();
