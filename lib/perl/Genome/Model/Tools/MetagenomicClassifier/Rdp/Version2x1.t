#!/usr/bin/env genome-perl

use strict;
use warnings;

use above "Genome";

use Test::More;

use_ok('Genome::Model::Tools::MetagenomicClassifier::Rdp::Version2x1') or die;

my $seq = {  
    id => 'S000002017 Pirellula staleyi',
    seq => 'AATGAACGTTGGCGGCATGGATTAGGCATGCAAGTCGTGCGCGATATGTAGCAATACATGGAGAGCGGCGAAAGGGAGAGTAATACGTAGGAACCTACCTTCGGGTCTGGGATAGCGGCGGGAAACTGCCGGTAATACCAGATGATGTTTCCGAACCAAAGGTGTGATTCCGCCTGAAGAGGGGCCTACGTCGTATTAGCTAGTTGGTAGGGTAATGGCCTACCAAGnCAAAGATGCGTATGGGGTGTGAGAGCATGCCCCCACTCACTGGGACTGAGACACTGCCCAGACACCTACGGGTGGCTGCAGTCGAGAATCTTCGGCAATGGGCGAAAGCCTGACCGAGCGATGCCGCGTGCGGGATGAAGGCCTTCGGGTTGTAAACCGCTGTCGTAGGGGATGAAGTGCTAGGGGGTTCTCCCTCTAGTTTGACTGAACCTAGGAGGAAGGnCCGnCTAATCTCGTGCCAGCAnCCGCGGTAATACGAGAGGCCCAnACGTTATTCGGATTTACTGGGCTTAAAGAGTTCGTAGGCGGTCTTGTAAGTGGGGTGTGAAATCCCTCGGCTCAACCGAGGAACTGCGCTCCAnACTACAAGACTTGAGGGGGATAGAGGTAAGCGGAACTGATGGTGGAGCGGTGAAATGCGTTGATATCATCAGGAACACCGGAGGCGAAGGCGGCTTACTGGGTCCTTTCTGACGCTGAGGAACGAAAGCTAGGGGAGCAnACGGGATTAGATACCCCGGTAGTCCTAnCCGTAAACGATGAGCACTGGACCGGAGCTCTGCACAGGGTTTCGGTCGTAGCGAAAGTGTTAAGTGCTCCGCCTGGGGAGTATGGTCGCAAGGCTGAAACTCAAAGGAATTGACGGGGGCTCACACAAGCGGTGGAGGATGTGGCTTAATTCGAGGCTACGCGAAGAACCTTATCCTAGTCTTGACATGCTTAGGAATCTTCCTGAAAGGGAGGAGTGCTCGCAAGAGAGCCTnTGCACAGGTGCTGCATGGCTGTCGTCAGCTCGTGTCGTGAGATGTCGGGTTAAGTCCCTTAACGAGCGAAACCCTnGTCCTTAGTTACCAGCGCGTCATGGCGGGGACTCTAAGGAGACTGCCGGTGTTAAACCGGAGGAAGGTGGGGATGACGTCAAGTCCTCATGGCCTTTATGATTAGGGCTGCACACGTCCTACAATnGTGCACACAAAGCGACGCAAnCTCGTGAGAGCCAGCTAAGTTCGGATTGCAGGCTGCAACTCGCCTGCATGAAGCTGGAATCGCTAGTAATCGCGGGTCAGCATACCGCGGTGAATGTGTTCCTGAGCCTTGTACACACCGCCCGTCAAGCCACGAAAGTGGGGGGGACCCAACAGCGCTGCCGTAACCGCAAGGAACAAGGCGCCTAAGGTCAACTCCGTGATTGGGACTAAGTCGTAACAAGGTAGCCGTAGGGGAACCTGCGGCTGGATCACCTCCTT',
};
$seq->{seq} = reverse $seq->{seq};
$seq->{seq} =~ tr/ATGC/TACG/;

my $classifier = Genome::Model::Tools::MetagenomicClassifier::Rdp::Version2x1->create(
    training_set => 4,# 4,6,broad
);
ok($classifier, 'Created rdp classifier');

my $classification = $classifier->classify($seq);
ok($classification, 'got classification from classifier');

my $i = 0;
my @taxa = (qw/ Root Bacteria Planctomycetes Planctomycetacia Planctomycetales Planctomycetaceae Pirellula /);
for my $rank (qw/ root domain phylum class order family genus /) {
    is($classification->{$rank}->{id}, $taxa[$i], 'taxon: '.$taxa[$i]);
    is($classification->{$rank}->{confidence}, '1.0', 'confidence: 1.0');
    $i++;
}
is($classification->{complemented}, 1, 'complemented');

# classify fails
eval{ $classifier->classify(); };
like($@, qr(No sequence given to classify), 'fail to classify w/ undef sequence');
diag($@);
eval{ $classifier->classify({seq => 'A'}); };
like($@, qr(Seq does not have an id:), 'fail to classify w/o sequence id');
ok(!$classifier->classify({ id => 'Short Seq', -seq => 'A' }), 'fail to classify short sequence');

# create fails
my $fail = Genome::Model::Tools::MetagenomicClassifier::Rdp::Version2x1->create();
ok(!$fail, 'fail w/o training set');
$fail = Genome::Model::Tools::MetagenomicClassifier::Rdp::Version2x1->create(training_set => 5);
ok(!$fail, 'fail w/ invalid training set');

done_testing();
