#!/usr/bin/env genome-perl
use above "Genome";
use Test::More tests => 5;
use Genome::Model::ClinSeq::TestData;
my $models = Genome::Model::ClinSeq::TestData::load();
my $m1id = 'b5d2ee5784e44a009886a0580f161b14';
my $m2id = 'de208ec52aee41819cc813a4287be955';

my $m1 = Genome::Model->get($m1id);
ok($m1, "got model 1");

my $m2 = Genome::Model->get($m2id);
ok($m2, "got model 2");

my $pair = Genome::Model::Pair->get(first => $m1, second => $m2);
ok($pair, "got a pair for models $m1id and $m2id");

is($pair->first, $m1, "first model is correct");
is($pair->second, $m2, "second model is correct");

