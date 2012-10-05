#!/usr/bin/env genome-perl

use strict;
use warnings;

use above "Genome";

require File::Compare;
use Test::More;

# Use
use_ok('Genome::Model::Tools::Sx::Trim::ByQual') or die;

# Create failures
ok(scalar(Genome::Model::Tools::Sx::Trim::ByQual->create->__errors__), 'Create w/o quality');
ok(scalar(Genome::Model::Tools::Sx::Trim::ByQual->create(quality => 'all')->__errors__), 'Create w/ quality => all');
ok(scalar(Genome::Model::Tools::Sx::Trim::ByQual->create(quality => 0)->__errors__), 'Create w/ quality => 0');

# Files
my $dir = $ENV{GENOME_TEST_INPUTS} . '/Genome-Model-Tools-Sx/TrimByQual';
my $in_fastq = $dir.'/in.fastq';
ok(-s $in_fastq, 'in fastq');
my $example_fastq = $dir.'/out.fastq';
ok(-s $example_fastq, 'example fastq');

my $tmp_dir = File::Temp::tempdir(CLEANUP => 1);
my $out_fastq = $tmp_dir.'/out.fastq';

# Ok
my $trimmer = Genome::Model::Tools::Sx::Trim::ByQual->create(
    input  => [ $in_fastq ],
    output => [ $out_fastq ],
    quality => 25,
);
ok($trimmer, 'create trimmer');
ok($trimmer->execute, 'execute trimmer');
is(File::Compare::compare($example_fastq, $out_fastq), 0, "fastq trimmed as expected");

print "gvimdiff $example_fastq $out_fastq\n"; <STDIN>;
done_testing();
exit;

