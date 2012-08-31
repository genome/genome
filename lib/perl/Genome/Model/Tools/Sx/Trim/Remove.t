#!/usr/bin/env genome-perl

use strict;
use warnings;

use above "Genome";

require File::Compare;
use Test::More;

# Use
use_ok('Genome::Model::Tools::Sx::Trim::Remove') or die;

# Create failures
ok(scalar(Genome::Model::Tools::Sx::Trim::Remove->create->__errors__), 'Create w/o length');
ok(scalar(Genome::Model::Tools::Sx::Trim::Remove->create(length => 'all')->__errors__), 'Create w/ length => all');
ok(scalar(Genome::Model::Tools::Sx::Trim::Remove->create(length => 0)->__errors__), 'Create w/ length => 0');

# Files
my $dir = $ENV{GENOME_TEST_INPUTS} . '/Genome-Model-Tools-Sx';
my $in_fastq = $dir.'/trimmer.input.fastq';
ok(-s $in_fastq, 'in fastq');
my $example_fastq = $dir.'/trimmer.remove.example.fastq';
ok(-s $example_fastq, 'example fastq');

my $tmp_dir = File::Temp::tempdir(CLEANUP => 1);
my $out_fastq = $tmp_dir.'/out.fastq';

# Ok
my $trimmer = Genome::Model::Tools::Sx::Trim::Remove->create(
    input  => [ $in_fastq ],
    output => [ $out_fastq ],
    length => 10,
);
ok($trimmer, 'create trimmer');
ok($trimmer->execute, 'execute trimmer');
is(File::Compare::compare($example_fastq, $out_fastq), 0, "fastq trimmed as expected");

#print "gvimdiff $out_fastq $example_fastq\n"; <STDIN>;
done_testing();
exit;

