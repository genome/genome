#!/usr/bin/env genome-perl

use strict;
use warnings;

use above 'Genome';

use Test::More skip_all => 'This test is incomplete';

my @cov    = qw[0 0 0 12 13 14 0 0 0 0 22 23 25 27 0 0 0];
my @seq    = qw[A G T  T  T  C G G A C  C  T  G  C C C C];
my $reflen = 17;

my $myRef   = Genome::Model::Tools::RefCov::Reference::GC->create();
$myRef->calculate_nucleotide_coverage(
    coverage => \@cov,
    sequence => \@seq,
);
print  Data::Dumper::Dumper($myRef);
$myRef->calculate_nucleotide_coverage(
    coverage => \@cov,
    sequence => \@seq,
);
print  Data::Dumper::Dumper($myRef);
exit;
print "GC REFLEN\t" . $myRef->GC_reflen_percent()   . "\n";
print "GC COV   \t" . $myRef->GC_covlen_percent()   . "\n";
print "GC UNCOV \t" . $myRef->GC_uncovlen_percent() . "\n";
print "\n";

print "AT REFLEN\t" . $myRef->AT_reflen_percent()   . "\n";
print "AT COV   \t" . $myRef->AT_covlen_percent()   . "\n";
print "AT UNCOV \t" . $myRef->AT_uncovlen_percent() . "\n";
print "\n";

print "G REFLEN\t" . $myRef->G_reflen_percent()   . "\n";
print "G COV   \t" . $myRef->G_covlen_percent()   . "\n";
print "G UNCOV \t" . $myRef->G_uncovlen_percent() . "\n";
print "\n";

print "N REFLEN\t" . $myRef->N_reflen_percent()   . "\n";
print "N COV   \t" . $myRef->N_covlen_percent()   . "\n";
print "N UNCOV \t" . $myRef->N_uncovlen_percent() . "\n";
print "\n";

for (1..17) { print $_ . "\t" }
print "\n";
print join ("\t", @cov) . "\n";
print join ("\t", @seq) . "\n";

# use YAML::Syck;
# print Dump( $myRef );
