#!/usr/bin/env genome-perl

use strict;
use warnings;

use above 'Genome';

require File::Compare;
use Test::More;

# use
use_ok('Genome::Model::Tools::Sx::Filter::ByMinLength') or die;

# create fail
ok(!Genome::Model::Tools::Sx::Filter::ByMinLength->execute(), 'execute w/o length');
ok(!Genome::Model::Tools::Sx::Filter::ByMinLength->execute(length => 'all'), 'execute w/ length => all');
ok(!Genome::Model::Tools::Sx::Filter::ByMinLength->execute(length => 0), 'execute w/ length => 0');

# Files
my $dir = $ENV{GENOME_TEST_INPUTS} . '/Genome-Model-Tools-Sx';
my $in_fastq = $dir.'/trimmer_bwa_style.example.fastq';
ok(-s $in_fastq, 'in fastq');
my $example_fastq = $dir.'/filter_by_length.example.fastq';
ok(-s $example_fastq, 'example fastq');

my $tmp_dir = File::Temp::tempdir(CLEANUP => 1);
my $out_fastq = $tmp_dir.'/out.fastq';

# Ok
my $filter = Genome::Model::Tools::Sx::Filter::ByMinLength->create(
    input  => [ $in_fastq ],
    output => [ $out_fastq ],
    length => 10,
);
ok($filter, 'create filter');
isa_ok($filter, 'Genome::Model::Tools::Sx::Filter::ByMinLength');
ok($filter->execute, 'execute filter');
is(File::Compare::compare($example_fastq, $out_fastq), 0, "fastq filtered as expected");

#print "$tmp_dir\n"; <STDIN>;
done_testing();
exit;

