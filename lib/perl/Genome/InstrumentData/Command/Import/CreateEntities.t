#!/usr/bin/env genome-perl

BEGIN {
    $ENV{UR_DBI_NO_COMMIT} = 1;
}

use strict;
use warnings;

use above 'Genome';

require File::Spec;
use Genome::Utility::Test;
use Test::Exception;
use Test::More;

my $class = 'Genome::InstrumentData::Command::Import::CreateEntities';
use_ok($class) or die;
my $data_dir = Genome::Utility::Test->data_dir_ok('Genome::InstrumentData::Command::Import', 'generate-cmds');
my $file = File::Spec->join($data_dir, 'create-entities.csv');

my %names_and_srs_sample_names = (
    Andromeda => 'SRS394779',
    Bwambale => 'SRS394801',
    Vincent => 'SRS394730',
);
for my $name ( keys %names_and_srs_sample_names ) {
    my $srs_name = $names_and_srs_sample_names{$name};
    my @libraries = Genome::Library->get(name => join('-', 'TGI', $name, $srs_name, 'extlibs'));
    is(@libraries, 0, "no $name libraries");
    my @samples = Genome::Sample->get(name => join('-', 'TGI', $name, $srs_name));
    is(@samples, 0, "no $name samples");
    my @individuals = Genome::Individual->get(name => join('-', 'TGI', $name));
    is(@individuals, 0, "no $name individuals");
}

# Fail w/o taxon existing
throws_ok(
    sub{ $class->execute(file => $file); },
    qr/Taxon does not exist for Pan troglodytes testerii\!/,
    'fails w/o taxon existing'
);
map { $_->delete } Genome::InstrumentData::Command::Import::Inputs->get;

# Create taxon
my $taxon = Genome::Taxon->create(name => 'Pan troglodytes testerii');
ok($taxon, 'create taxon');

my $cmd = $class->execute(
    file => $file,
);
ok($cmd->result, 'execute');
for my $name ( keys %names_and_srs_sample_names ) {
    my $srs_name = $names_and_srs_sample_names{$name};
    my @libraries = Genome::Library->get(name => join('-', 'TGI', $name, $srs_name, 'extlibs'));
    is(@libraries, 1, "one $name libraries");
    my @samples = Genome::Sample->get(name => join('-', 'TGI', $name, $srs_name));
    is(@samples, 1, "one $name samples");
    is($libraries[0]->sample, $samples[0], 'correct sample library association');
    my @individuals = Genome::Individual->get(name => join('-', 'TGI', $name));
    is(@individuals, 1, "one $name individuals");
    is($samples[0]->source, $individuals[0], 'correct individual sample association');
}
map { $_->delete } Genome::InstrumentData::Command::Import::Inputs->get;

# Rexecute to make sure things are not recreated
$cmd = $class->execute(
    file => $file,
);
ok($cmd->result, 'execute');
for my $name ( keys %names_and_srs_sample_names ) {
    my $srs_name = $names_and_srs_sample_names{$name};
    my @libraries = Genome::Library->get(name => join('-', 'TGI', $name, $srs_name, 'extlibs'));
    is(@libraries, 1, "one $name libraries");
    my @samples = Genome::Sample->get(name => join('-', 'TGI', $name, $srs_name));
    is(@samples, 1, "one $name samples");
    is($libraries[0]->sample, $samples[0], 'correct sample library association');
    my @individuals = Genome::Individual->get(name => join('-', 'TGI', $name));
    is(@individuals, 1, "one $name individuals");
    is($samples[0]->source, $individuals[0], 'correct individual sample association');
}
map { $_->delete } Genome::InstrumentData::Command::Import::Inputs->get;

# Fails
throws_ok(
    sub{ $class->execute(file => File::Spec->join($data_dir, 'missing-required-property.csv')); },
    qr/No taxon for/,
    'failed when missing taxon',
);

done_testing();
