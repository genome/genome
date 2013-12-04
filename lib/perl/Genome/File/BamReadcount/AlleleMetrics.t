#!/usr/bin/env perl

use above 'Genome';
use Test::More;

use strict;
use warnings;

use_ok("Genome::File::BamReadcount::AlleleMetrics");


subtest 'Normal allele metrics parsing' => sub {
    my $test_allele_metrics = "T:96:36.73:26.54:33.20:50:46:0.49:0.02:17.53:92:0.41:229.41:0.40"; 
    my $metric = Genome::File::BamReadcount::AlleleMetrics->new($test_allele_metrics);
    ok($metric, "Metric created without error");

    is($metric->allele, 'T', "Allele is correct");
    is($metric->count, 96, "Count is correct");
    is($metric->avg_mapq, 36.73, "Average Mapping Quality is correct");
    is($metric->avg_bq, 26.54, "Average Base Quality is correct");
    is($metric->avg_se_mapq, '33.20', "Average Single-ended Mapping Quality is correct");
    is($metric->num_plus_strand, 50, "Number of reads on plus strand is correct");
    is($metric->num_minus_strand, 46, "Number of reads on minus strand is correct");
    is($metric->avg_pos_as_fraction, 0.49, "Average position as a fraction of readlength is correct");
    is($metric->avg_num_mismatches_as_fraction, 0.02, "Average number of mismatches as a fraction of readlength is correct");
    is($metric->avg_sum_mismatch_qualities, 17.53, "Average mismatch quality sum is correct");
    is($metric->num_q2_containing_reads, 92, "Number of Q2 containing reads is correct");
    is($metric->avg_distance_to_q2_start_in_q2_reads, 0.41, "Average distance to Q2 start is correct");
    is($metric->avg_clipped_length, 229.41, "Average clipped readlength is correct");
    is($metric->avg_distance_to_effective_3p_end, '0.40', "Average distrance to effective 3' end is correct");

    ok(!$metric->is_indel, "Not an indel");
    ok(!$metric->is_del, "Not a deletion");
    ok(!$metric->is_ins, "Not an insertion");
};

subtest "Indel metric parsing" => sub {
    my $deletion_metrics = "-GCT:96:36.73:26.54:0.00:50:46:0.49:0.02:17.53:92:0.41:229.41:0.40"; 
    my $insertion_metrics = "+G:96:36.73:26.54:0.00:50:46:0.49:0.02:17.53:92:0.41:229.41:0.40"; 

    my $deletion_metric = Genome::File::BamReadcount::AlleleMetrics->new($deletion_metrics);
    ok($deletion_metric, "Deletion metric created without error");
    ok($deletion_metric->is_indel, "Is an indel");
    ok($deletion_metric->is_del, "Is a deletion");
    ok(!$deletion_metric->is_ins, "Is not an insertion");
    is($deletion_metric->allele, "-GCT", "Deletion allele matches");

    my $insertion_metric = Genome::File::BamReadcount::AlleleMetrics->new($insertion_metrics);
    ok($insertion_metric, "Insertion metric created without error");
    ok($insertion_metric->is_indel, "Is an indel");
    ok(!$insertion_metric->is_del, "Is not an insertion");
    ok($insertion_metric->is_ins, "Is an insertion");
    is($insertion_metric->allele, "+G", "Insertion allele matches");
};

subtest "Malformed metric" => sub {
    my $truncated_allele_metrics = "T:96:36.73:26.54:"; 
    eval {
        my $metric = Genome::File::BamReadcount::AlleleMetrics->new($truncated_allele_metrics);
    };
    ok($@, "Truncated metric throws an error");
};
done_testing();
