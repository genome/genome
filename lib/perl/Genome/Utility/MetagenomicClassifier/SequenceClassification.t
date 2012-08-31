#!/usr/bin/env genome-perl

use strict;
use warnings;

use above 'Genome';

use Bio::Taxon;
use Test::More;

use_ok('Genome::Utility::MetagenomicClassifier::SequenceClassification') or die;
use_ok('Genome::Utility::MetagenomicClassifier');

my @taxa;
my $string = 'Root:1.0;Bacteria:1.0;Eubacteria:1.0;Bacteroidetes:0.99;Bacteroidetes:0.82;Bacteroidales:0.82;Rikenellaceae:0.78;Alistipes:0.68;Alistipes carmichaelli:0.68';
my @ranks = Genome::Utility::MetagenomicClassifier->taxonomic_ranks;
unshift @ranks, 'root';
for my $assignment ( split(';', $string) ) {
    my ($name, $conf) = split(':', $assignment);
    push @taxa, Genome::Utility::MetagenomicClassifier->create_taxon(
        id => $name,
        rank => shift(@ranks),
        tags => {
            confidence => $conf,
        },
        ancestor => ( @taxa ? $taxa[$#taxa] : undef ),
    );
}

my %params = (
    name => 'U_PR-aab10d09',
    complemented => 0,
    classifier => 'rdp',
    taxon => $taxa[0],
);
my $sequence_classification = Genome::Utility::MetagenomicClassifier::SequenceClassification->new(%params);
ok($sequence_classification, 'create sequence classification');

for my $key ( keys %params ) {
    my $method = 'get_'.$key;
    can_ok($sequence_classification, $method);
    is_deeply([$sequence_classification->$method], [$params{$key}], $key);
}
can_ok($sequence_classification, 'get_taxa');

for my $rank ( 'root', Genome::Utility::MetagenomicClassifier->taxonomic_ranks ) {
    # taxon
    my $get_taxon_method = 'get_'.$rank.'_taxon';
    can_ok($sequence_classification, $get_taxon_method);
    my $taxon = $sequence_classification->$get_taxon_method;
    ok($taxon, "Got $rank taxon");
    is($taxon->rank, $rank, "Taxon is $rank");
    # name
    my $get_name_method = 'get_'.$rank;
    my $name = $sequence_classification->$get_name_method;
    ok($name, "Got $rank name ($name) for taxon");
    is($taxon->id, $name, "Taxon name and $get_name_method match");
    # name and confidence
    my $get_conf_method = 'get_'.$rank.'_confidence';
    my $conf = $sequence_classification->$get_conf_method;
    ok($conf, "Got confidence ($conf) for $rank with $get_conf_method");
    my ($conf_from_taxon) = $taxon->get_tag_values('confidence');
    is($conf, $conf_from_taxon, "Confidence from $get_conf_method and taxon match");
}

# Check that these private methods do not return stuff that doesn't exist
ok(!$sequence_classification->_get_taxon_for_rank('blah'), "As expected - no blah taxon");
ok(!$sequence_classification->_get_taxon_name_for_rank('blah'), "As expected - no blah taxon name");
ok(!$sequence_classification->_get_taxon_confidence_for_rank('blah'), "As expected - no blah confidence");

done_testing();
exit;

