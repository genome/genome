#!/usr/bin/env genome-perl

use strict;
use warnings;

use above 'Genome';

require File::Compare;
use Test::More;

# use
use_ok('Genome::Model::Tools::Sx::Filter::ByMaxNs') or die;

# create fail
ok(!Genome::Model::Tools::Sx::Filter::ByMaxNs->execute(), 'execute w/o maximum');
ok(!Genome::Model::Tools::Sx::Filter::ByMaxNs->execute(maximum => 'all'), 'execute w/ maximum => all');
ok(!Genome::Model::Tools::Sx::Filter::ByMaxNs->execute(maximum => -1), 'execute w/ maximum less than 0');

# Files
my $dir = $ENV{GENOME_TEST_INPUTS} . '/Genome-Model-Tools-Sx/FilterByMaxNs';
my $in_fwd_fastq = $dir.'/in.fwd.fastq';
ok(-s $in_fwd_fastq, 'in fastq');
my $in_rev_fastq = $dir.'/in.rev.fastq';
ok(-s $in_rev_fastq, 'in fastq');
my $expected_out_fwd_fastq = $dir.'/out.fwd.fastq';
ok(-s $expected_out_fwd_fastq, 'expected fwd fastq');
my $expected_out_rev_fastq = $dir.'/out.rev.fastq';
ok(-s $expected_out_rev_fastq, 'expected rev fastq');

my $tmp_dir = File::Temp::tempdir(CLEANUP => 1);
my $out_fwd_fastq = $tmp_dir.'/out.fwd.fastq';
my $out_rev_fastq = $tmp_dir.'/out.rev.fastq';

# Ok
my $filter = Genome::Model::Tools::Sx::Filter::ByMaxNs->create(
    input  => [ $in_fwd_fastq, $in_rev_fastq ],
    output => [ $out_fwd_fastq.':name=fwd', $out_rev_fastq.':name=rev' ],
    maximum => 15,
);
ok($filter, 'create filter');
isa_ok($filter, 'Genome::Model::Tools::Sx::Filter::ByMaxNs');
ok($filter->execute, 'execute filter');
is(File::Compare::compare($expected_out_fwd_fastq, $out_fwd_fastq), 0, "fwd fastq filtered as expected");
is(File::Compare::compare($expected_out_rev_fastq, $out_rev_fastq), 0, "rev fastq filtered as expected");

#print "$tmp_dir\n"; <STDIN>;
done_testing();
