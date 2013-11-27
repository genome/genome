#!/usr/bin/env perl

use above 'Genome';
use Test::More;

use strict;
use warnings;

use_ok("Genome::File::BamReadcount::Entry");

subtest 'Regular reporting' => sub {
    my $test_string = "21	10400136	T	98	=:0:0.00:0.00:0.00:0:0:0.00:0.00:0.00:0:0.00:0.00:0.00	A:1:42.00:2.00:0.00:0:1:0.78:0.04:14.00:1:0.13:201.00:0.13	C:1:30.00:2.00:0.00:1:0:0.55:0.13:42.00:0:nan:199.00:0.22	G:0:0.00:0.00:0.00:0:0:0.00:0.00:0.00:0:0.00:0.00:0.00	T:96:36.73:26.54:0.00:50:46:0.49:0.02:17.53:92:0.41:229.41:0.40	N:0:0.00:0.00:0.00:0:0:0.00:0.00:0.00:0:0.00:0.00:0.00";
    my $new_entry = Genome::File::BamReadcount::Entry->new($test_string);
    ok($new_entry, "Created regular output entry");
    is($new_entry->chromosome, "21", "Chromosome correct");
    is($new_entry->position, "10400136", "Position correct");
    is($new_entry->ref_base, "T", "Reference base correct");
    is($new_entry->depth, 98, "Depth correct");
    is($new_entry->libraries, undef, "Libraries undefined");
    is($new_entry->num_libraries, 0, "Number of libraries is 0");
    ok(!$new_entry->has_per_library, "Does not report as having per library metrics");
    my $allele_metrics = $new_entry->metrics_for("T");
    ok($allele_metrics, "Able to retrieve metrics for base 'T'");
    is($allele_metrics->allele, 'T', "Reported allele matches expected");
};

subtest 'Per-lib reporting' => sub {
    my $per_lib_test_string = "21	10402985	G	344	Solexa-135852	{	=:0:0.00:0.00:0.00:0:0:0.00:0.00:0.00:0:0.00:0.00:0.00	A:2:57.50:16.00:0.00:2:0:0.35:0.02:42.50:2:0.51:231.50:0.51	C:0:0.00:0.00:0.00:0:0:0.00:0.00:0.00:0:0.00:0.00:0.00	G:155:53.61:25.61:0.00:90:65:0.51:0.01:11.10:134:0.39:236.59:0.38	T:0:0.00:0.00:0.00:0:0:0.00:0.00:0.00:0:0.00:0.00:0.00	N:0:0.00:0.00:0.00:0:0:0.00:0.00:0.00:0:0.00:0.00:0.00	+A:20:52.55:0.00:0.00:12:8:0.53:0.02:25.00:18:0.41:237.70:0.39	}	Solexa-135853	{	=:0:0.00:0.00:0.00:0:0:0.00:0.00:0.00:0:0.00:0.00:0.00	A:1:60.00:2.00:0.00:0:1:0.76:0.10:46.00:0:-nan:250.00:0.62	C:0:0.00:0.00:0.00:0:0:0.00:0.00:0.00:0:0.00:0.00:0.00	G:186:52.12:26.33:0.29:99:87:0.50:0.01:12.29:161:0.39:238.74:0.40	T:0:0.00:0.00:0.00:0:0:0.00:0.00:0.00:0:0.00:0.00:0.00	N:0:0.00:0.00:0.00:0:0:0.00:0.00:0.00:0:0.00:0.00:0.00	}";

    my $second_entry = Genome::File::BamReadcount::Entry->new($per_lib_test_string);
    ok($second_entry, "Created per-lib output entry");
    is($second_entry->chromosome, "21", "Chromosome correct");
    is($second_entry->position, 10402985, "Position correct");
    is($second_entry->ref_base, "G", "Reference base correct");
    is($second_entry->depth, 344, "Depth correct");
    is($second_entry->num_libraries, 2, "Num libraries correct");
    ok($second_entry->has_per_library, "Reports as having per-library");
    my @libmetrics = $second_entry->libraries;
    my @libnames = map {$_->name} @libmetrics;
    is_deeply(\@libnames, [qw(Solexa-135852 Solexa-135853)], "Names report correctly");
};

subtest 'Truncated entry' => sub {
    my $bad_entry_string = "X	20000	C";
    eval {
        my $entry = Genome::File::BamReadcount::Entry->new($bad_entry_string);
    };
    ok($@, "Entry throws on truncated input");
};

done_testing();
