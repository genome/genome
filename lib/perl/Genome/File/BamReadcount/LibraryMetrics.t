#!/usr/bin/env perl

use above 'Genome';
use Test::More;

use strict;
use warnings;

use_ok("Genome::File::BamReadcount::LibraryMetrics");

subtest 'Normal library metrics parsing' => sub {
    my $lib_test_string = "Whatever	{	=:0:0.00:0.00:0.00:0:0:0.00:0.00:0.00:0:0.00:0.00:0.00	A:2:57.50:16.00:0.00:2:0:0.35:0.02:42.50:2:0.51:231.50:0.51	C:0:0.00:0.00:0.00:0:0:0.00:0.00:0.00:0:0.00:0.00:0.00	G:155:53.61:25.61:0.00:90:65:0.51:0.01:11.10:134:0.39:236.59:0.38	T:0:0.00:0.00:0.00:0:0:0.00:0.00:0.00:0:0.00:0.00:0.00	N:0:0.00:0.00:0.00:0:0:0.00:0.00:0.00:0:0.00:0.00:0.00	+A:20:52.55:0.00:0.00:12:8:0.53:0.02:25.00:18:0.41:237.70:0.39	}";
    my $libmetric = Genome::File::BamReadcount::LibraryMetrics->new($lib_test_string);
    ok($libmetric, "LibraryMetric created without error");
    is($libmetric->name, "Whatever", "Library name is correct");
    is($libmetric->depth, 177, "Library depth is correct");
    my @alleles = $libmetric->alleles;
    is_deeply(\@alleles, [qw( = A C G T N +A )], "Alleles returned as expected");
    
    subtest 'Test AlleleMetric interface' => sub {
        my $allele_metric = $libmetric->metrics_for('+A');
        ok($allele_metric, "Retrieved allele metric");
        eval {
            my $bad_allele_metric = $libmetric->metrics_for('Prepare to die, obviously');
        };
        ok($@, "Exception thrown when trying to retrieve an allele that doesn't exist");
        is($allele_metric->allele, '+A', "Allele as expected");
        is($allele_metric->count, 20, "Allele count as expected");
    };
};

subtest 'Truncated library metrics parsing' => sub {
    my $truncated_test_string = "Cool	{	=:0:0.00:0.00:0.00:0:0:0.00:0.00:0.00:0:0.00:0.00:0.00	A:2:57.50:16.00:0.00:2:0:0.35:0.02:42.50:2:0.51:231.50:0.51	C:0:0.00:0.00:0.00:0:0:0.00:0.00:0.00:0:0.00:0.00:0.00	G:155:53.61:25.61:0.00:90:65:0.51:0.01:11.10:134:0.39:236.59:0.38";
    eval {
        my $libmetric = Genome::File::BamReadcount::LibraryMetric->new($truncated_test_string);
    };
    ok($@, "Truncated metric throws an error");
};

done_testing();
