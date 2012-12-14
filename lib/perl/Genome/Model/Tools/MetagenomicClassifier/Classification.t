#!/usr/bin/env genome-perl

use strict;
use warnings;

use above 'Genome';

use Data::Dumper 'Dumper';
use File::Compare 'compare';
require Storable;
require File::Temp;
use Test::More;

use_ok('Genome::Model::Tools::MetagenomicClassifier::ClassificationReader') or die;
use_ok('Genome::Model::Tools::MetagenomicClassifier::ClassificationWriter') or die;
use_ok('Genome::Model::Tools::MetagenomicClassifier::ClassificationComposition') or die;

my $base_dir = $ENV{GENOME_TEST_INPUTS} . '/Genome-Utility-MetagenomicClassifier-SequenceClassification-Writer/';
my $tmp_dir = File::Temp::tempdir(CLEANUP => 1);

# read/write hmp_fix_ranks
my $reader = Genome::Model::Tools::MetagenomicClassifier::ClassificationReader->create(
    file => $base_dir.'/classifications.hmp_fix_ranks',
);
ok($reader, 'create reader for hmp_fix_ranks');
my $writer = Genome::Model::Tools::MetagenomicClassifier::ClassificationWriter->create(
    file => $tmp_dir.'/classifications.hmp_fix_ranks',
    format => 'hmp_fix_ranks',
);
ok($writer, 'create writer for hmp_fix_ranks');
my $comp = Genome::Model::Tools::MetagenomicClassifier::ClassificationComposition->create(
    confidence_threshold => 0.8,
);
ok($comp, 'create classification composition');

my $classifications = 0;
while ( my $classification = $reader->read ) {
    $classifications++;
    ok($writer->write($classification), 'write classification');
    ok($comp->add_classification($classification), 'add classification');
}
is($classifications, 6, 'read/write 6 classifications');
is(File::Compare::compare($reader->file, $writer->file), 0, 'files match for hmp_fix_ranks');

my $confident_classifications = $comp->confident_classifications;
is(@$confident_classifications, 6, 'confident classifications');
my %counts = $comp->get_counts_for_domain_down_to_rank('bacteria', 'genus');
is_deeply(\%counts, _expected_counts(), 'counts for composition');

# read/write hmp_all_ranks
$reader = Genome::Model::Tools::MetagenomicClassifier::ClassificationReader->create(
    file => $base_dir.'/classifications.hmp_all_ranks',
);
ok($reader, 'create reader for hmp_fix_ranks');
$writer = Genome::Model::Tools::MetagenomicClassifier::ClassificationWriter->create(
    file => $tmp_dir.'/classifications.hmp_all_ranks',
    format => 'hmp_all_ranks',
);
ok($writer, 'create writer for hmp_all_ranks');

$classifications = 0;
while ( my $classification = $reader->read ) {
    $classifications++;
    $writer->write($classification);
}
is($classifications, 6, 'read/write 6 classifications');
is(File::Compare::compare($reader->file, $writer->file), 0, 'files match for hmp_fix_ranks');

#print 'gvimdiff '.$reader->file.' '.$writer->file."\n"; <STDIN>;
done_testing();


sub _expected_counts {
    return {
        'Bacteria:Firmicutes:Clostridia:Clostridiales:Veillonellaceae:Dialister' => {
            'domain' => 1,
            'order' => 1,
            'genus' => 1,
            'class' => 1,
            'phylum' => 1,
            'total' => 1,
            'family' => 1
        },
        'Bacteria:Firmicutes:Clostridia:Clostridiales:Ruminococcaceae:Sporobacter' => {
            'domain' => 1,
            'order' => 1,
            'genus' => 0,
            'class' => 1,
            'phylum' => 1,
            'total' => 1,
            'family' => 0
        },
        'Bacteria:Verrucomicrobia:Verrucomicrobia_no_class:Verrucomicrobia_no_order:Verrucomicrobia_no_family:Akkermansia' => {
            'domain' => 1,
            'order' => 1,
            'genus' => 1,
            'class' => 1,
            'phylum' => 1,
            'total' => 1,
            'family' => 1
        },
        'Bacteria:Firmicutes:Clostridia:Clostridiales:Ruminococcaceae:Papillibacter' => {
            'domain' => 1,
            'order' => 1,
            'genus' => 0,
            'class' => 1,
            'phylum' => 1,
            'total' => 1,
            'family' => 1
        },
        'Bacteria:Firmicutes:Clostridia:Clostridiales:Clostridiaceae:Oxobacter' => {
            'domain' => 1,
            'order' => 0,
            'genus' => 0,
            'class' => 0,
            'phylum' => 0,
            'total' => 1,
            'family' => 0
        },
        'Bacteria:Bacteroidetes:Bacteroidetes:Bacteroidales:Rikenellaceae:Alistipes' => {
            'domain' => 1,
            'order' => 1,
            'genus' => 0,
            'class' => 1,
            'phylum' => 1,
            'total' => 1,
            'family' => 0
        }
    };
}
