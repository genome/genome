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
        # these are deprecated
        name_regexp => { default_value => '(TGI\-[\w\d]+)\-[\w\d]+', },
        taxon_name => { default_value => 'human', },
    },
};

my $sample_name = 'TGI-INDINDIVIDUAL-SAMPLE';
my $sample = Genome::Sample->__define__(
    name => $sample_name,
);
ok($sample, 'defined sample');

# Fail w/ attrs missing value
my $fail = SampleImporter->execute(
    name => $sample_name,
    nomenclature => 'TGI',
    individual_attributes => [qw/ awesomeness= /],
);
ok(!$fail->result, 'execute sample importer failed w/o taxon');
is($fail->error_message, 'Individual attribute label (awesomeness) does not have a value!', 'correct error message');

$fail = SampleImporter->execute(
    name => $sample_name,
    nomenclature => 'TGI',
    sample_attributes => [qw/ awesomeness= /],
);
ok(!$fail->result, 'execute sample importer failed w/o taxon');
is($fail->error_message, 'Sample attribute label (awesomeness) does not have a value!', 'correct error message');

$fail = SampleImporter->execute(
    name => $sample_name,
    nomenclature => 'TGI',
    library_attributes => [qw/ awesomeness= /],
);
ok(!$fail->result, 'execute sample importer failed w/o taxon');
is($fail->error_message, 'Library attribute label (awesomeness) does not have a value!', 'correct error message');

# Success
my $cmd1 = SampleImporter->execute(name => $sample_name, nomenclature => 'TGI');
ok($cmd1->result, 'execute sample importer for existing sample');
ok($cmd1->_individual, 'set individual');
is($cmd1->_sample, $sample, 'set existing sample');
ok($cmd1->_library, 'set library');

# Library already imported = success
my $cmd2 = SampleImporter->execute(name => $sample_name, nomenclature => 'TGI'); # taxon not needed since ind exists
ok($cmd2->result, 'execute sample importer for existing sample');
ok(!$cmd2->_individual, 'did not set individual');
is($cmd2->_sample, $sample, 'found existing sample');
is($cmd2->_library, $cmd1->_library, 'found existing library');

done_testing();
