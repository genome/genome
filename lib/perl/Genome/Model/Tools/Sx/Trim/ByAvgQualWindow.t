#!/usr/bin/env genome-perl

use strict;
use warnings;

use above "Genome";

require File::Compare;
use Test::More;

# Use
use_ok('Genome::Model::Tools::Sx::Trim::ByAvgQualWindow') or die;

# Create failures
ok(scalar(Genome::Model::Tools::Sx::Trim::ByAvgQualWindow->create->__errors__), 'Create w/o quality and window');
ok(scalar(Genome::Model::Tools::Sx::Trim::ByAvgQualWindow->create(quality => 'all', window => 50)->__errors__), 'Create w/ invalid quality (all)');
ok(scalar(Genome::Model::Tools::Sx::Trim::ByAvgQualWindow->create(quality => 0, window => 50)->__errors__), 'Create w/ invalid quality (0)');
ok(scalar(Genome::Model::Tools::Sx::Trim::ByAvgQualWindow->create(quality => 25, window => 'all')->__errors__), 'Create w/ invalid window (all)');
ok(scalar(Genome::Model::Tools::Sx::Trim::ByAvgQualWindow->create(quality => 25, window => 0)->__errors__), 'Create w/ invalid window (0)');

# Files
my $dir = $ENV{GENOME_TEST_INPUTS} . '/Genome-Model-Tools-Sx';
my $in_fastq = $dir.'/trimmer.in.fastq';
ok(-s $in_fastq, 'in fastq');
my $example_fastq = $dir.'/trimmer_by_avg_qual_window.example2.fastq';
ok(-s $example_fastq, 'example fastq');

my $tmp_dir = File::Temp::tempdir(CLEANUP => 1);
my $out_fastq = $tmp_dir.'/out.fastq';

# Ok
my $trimmer = Genome::Model::Tools::Sx::Trim::ByAvgQualWindow->create(
    input  => [ $in_fastq ],
    output => [ $out_fastq ],
    quality => 25,
    window => 50,
);
ok($trimmer, 'create trimmer');
ok($trimmer->execute, 'execute trimmer');
is(File::Compare::compare($example_fastq, $out_fastq), 0, "fastq trimmed as expected");

#print "gvimdiff $in_fastq $out_fastq\n"; print "gvimdiff $example_fastq $out_fastq\n"; <STDIN>;
done_testing();
