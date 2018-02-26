#!/usr/bin/env genome-perl

use strict;
use warnings;

use above 'Genome';

use Data::Dumper 'Dumper';
use File::Compare 'compare';
require Storable;
require File::Temp;
use Test::More;

use_ok('Genome::Utility::MetagenomicClassifier::SequenceClassification::Writer');

# Seq classifications to write
use_ok('Genome::Utility::MetagenomicClassifier::SequenceClassification');
my $base_dir = Genome::Config::get('test_inputs') . '/Genome-Utility-MetagenomicClassifier-SequenceClassification-Writer/';
my $classifications_stor_file = $base_dir.'/classifications.stor';
my $classifications = Storable::retrieve($classifications_stor_file);
ok($classifications, 'Loaded classifications');

# Good w/ hmp_fix_ranks (default)
my $tmpdir = File::Temp::tempdir(CLEANUP => 1);
my $writer = Genome::Utility::MetagenomicClassifier::SequenceClassification::Writer->create(
    output => $tmpdir.'/classification.fix_ranks',
);
ok($writer, 'Created sequence classification writer for default fix ranks');
for my $classification ( @$classifications ) {
    $writer->write_one($classification);
}
my $hmp_fix_ranks_file = $base_dir.'/classifications.hmp_fix_ranks';
ok(-s $hmp_fix_ranks_file, 'HMP fix ranks file exists');
is(
    compare($writer->get_original_output, $hmp_fix_ranks_file), 
    0,
    'Generated and expected classification files match for fix ranks',
);
#print "diff $hmp_fix_ranks_file ".$writer->get_original_output."\n";<STDIN>;

# Good w/ hmp_all_ranks
$writer = Genome::Utility::MetagenomicClassifier::SequenceClassification::Writer->create(
    output => $tmpdir.'/classification.all_ranks',
    format => 'hmp_all_ranks',
);
ok($writer, 'Created sequence classification writer for all ranks');
for my $classification ( @$classifications ) {
    $writer->write_one($classification);
}
my $hmp_all_ranks_file = $base_dir.'/classifications.hmp_all_ranks';
ok(-s $hmp_all_ranks_file, 'HMP all ranks file exists');
is(
    compare($writer->get_original_output, $hmp_all_ranks_file), 
    0,
    'Generated and expected classification files match for all ranks',
);
#print "diff $hmp_all_ranks_file ".$writer->get_original_output."\n";<STDIN>;

# Fail
$writer = Genome::Utility::MetagenomicClassifier::SequenceClassification::Writer->create(
    output => $tmpdir.'/some_file.out',
    format => 'h',
    
);
ok(!$writer, 'Failed as expected - create w/ invalid format');

done_testing();
