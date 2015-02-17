#!/usr/bin/env genome-perl

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
    },
};

my $taxon = Genome::Taxon->__define__(name => 'almost human');
ok($taxon, 'defined taxon');
my $sample_name = 'TGI-INDINDIVIDUAL-SAMPLE';

# Fail w/o taxon
throws_ok(
    sub{ SampleImporter->execute(name => $sample_name, nomenclature => 'TGI'); },
    qr/No taxon given and an individual must be created\. Please specify one\./,
    'execute sample importer failed w/o taxon',
);

done_testing();
