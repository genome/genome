#! /gsc/bin/perl

use strict;
use warnings;

use above 'Genome';

use Test::More;

use_ok('Genome::Site::TGI::Synchronize::Classes::Dictionary') or die;

my $class_mapping = Genome::Site::TGI::Synchronize::Classes::Dictionary->get;
ok($class_mapping, 'get class mapping');

my @entities_to_sync = $class_mapping->entity_names;
is(@entities_to_sync, 10, 'entity names');

ok(!eval{$class_mapping->lims_class_for_entity_name;}, 'failed to get lims class for undef entity');
is($class_mapping->error_message, 'No entity name given to get LIMS class!', 'correct error');
ok(!$class_mapping->lims_class_for_entity_name('unknown'), 'no lims class for "unknown"');
is($class_mapping->error_message, 'No LIMS class for entity! unknown', 'correct error');
my $lims_class_for_taxon = $class_mapping->lims_class_for_entity_name('taxon');
is($lims_class_for_taxon, 'Genome::Site::TGI::Synchronize::Classes::OrganismTaxon', 'correct lims class for taxon');

ok(!eval{$class_mapping->genome_class_for_entity_name;}, 'failed to get genome class for undef entity');
is($class_mapping->error_message, 'No entity name given to get LIMS class!', 'correct error');
ok(!$class_mapping->genome_class_for_entity_name('unknown'), 'no genome class for "unknown"');
is($class_mapping->error_message, 'No LIMS class for entity! unknown', 'correct error');
my $genome_class_for_taxon = $class_mapping->genome_class_for_entity_name('taxon');
is($genome_class_for_taxon, 'Genome::Taxon', 'correct genome class for taxon');

done_testing();
