#!/usr/bin/env genome-perl

use strict;
use warnings;

use above 'Genome';

require File::Compare;
use Test::More;

# Use
use_ok('Genome::Model::Tools::Sx::Trim::BwaStyle') or die;

# Create fails
ok(
    !Genome::Model::Tools::Sx::Trim::BwaStyle->execute(trim_qual_level => 'pp'),
    'create w/ trim_qual_level => pp'
);
ok(
    !Genome::Model::Tools::Sx::Trim::BwaStyle->execute(trim_qual_level => -1),
    'create w/ trim_qual_level => -1'
);

# Files
my $dir = $ENV{GENOME_TEST_INPUTS} . '/Genome-Model-Tools-Sx';
my $in_fastq = $dir.'/trimmer.in.fastq';
ok(-s $in_fastq, 'in fastq');
my $example_fastq = $dir.'/trimmer_bwa_style.example.fastq';
ok(-s $example_fastq, 'example fastq');

my $tmp_dir = File::Temp::tempdir(CLEANUP => 1);
my $out_fastq = $tmp_dir.'/out.fastq';

# Ok
my $trimmer = Genome::Model::Tools::Sx::Trim::BwaStyle->create(
    input  => [ $in_fastq ],
    output => [ $out_fastq ],
    trim_qual_level => 10,
);
ok($trimmer, 'create trimmer');
isa_ok($trimmer, 'Genome::Model::Tools::Sx::Trim::BwaStyle');
ok($trimmer->execute, 'execute trimmer');
is(File::Compare::compare($example_fastq, $out_fastq), 0, "fastq trimmed as expected");

#print "diff $example_fastq $out_fastq\n"; <STDIN>;
done_testing();
