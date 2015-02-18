#!/usr/bin/env genome-perl

BEGIN {
    $ENV{UR_DBI_NO_COMMIT} = 1;
    $ENV{UR_COMMAND_DUMP_STATUS_MESSAGES} = 1;
}

use strict;
use warnings;

use above 'Genome';
use Test::Exception;
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

# Fail w/o taxon
throws_ok(
    sub{ SampleImporter->execute(name => $sample_name, nomenclature => 'TGI'); },
    qr/No taxon given and an individual must be created\. Please specify one\./,
    'execute sample importer failed w/o taxon',
);
my $sample = Genome::Sample->__define__(
    name => $sample_name,
);
ok($sample, 'defined sample');

done_testing();
