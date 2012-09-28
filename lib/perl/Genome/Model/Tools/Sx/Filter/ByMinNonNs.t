#!/usr/bin/env genome-perl

use strict;
use warnings;

use above 'Genome';

require File::Compare;
use Test::More;

# use
use_ok('Genome::Model::Tools::Sx::Filter::ByMinNonNs') or die;

# create fail
ok(!Genome::Model::Tools::Sx::Filter::ByMinNonNs->execute(), 'execute w/o minimum');
ok(!Genome::Model::Tools::Sx::Filter::ByMinNonNs->execute(minimum => 'all'), 'execute w/ minimum => all');
ok(!Genome::Model::Tools::Sx::Filter::ByMinNonNs->execute(minimum => -1), 'execute w/ minimum less than 0');

# Files
my $dir = $ENV{GENOME_TEST_INPUTS} . '/Genome-Model-Tools-Sx/FilterByMinNonNs';
my $in_fwd_fastq = $dir.'/in.fwd.fastq';
ok(-s $in_fwd_fastq, 'in fastq');
my $in_rev_fastq = $dir.'/in.rev.fastq';
ok(-s $in_rev_fastq, 'in fastq');
my $expected_out_fastq = $dir.'/out.fastq';
ok(-s $expected_out_fastq, 'expected fastq');

my $tmp_dir = File::Temp::tempdir(CLEANUP => 1);
my $out_fastq = $tmp_dir.'/out.rev.fastq';

# Ok
my $filter = Genome::Model::Tools::Sx::Filter::ByMinNonNs->create(
    input  => [ $in_fwd_fastq, $in_rev_fastq ],
    output => [ $out_fastq ],
    minimum => 60,
);
ok($filter, 'create filter');
isa_ok($filter, 'Genome::Model::Tools::Sx::Filter::ByMinNonNs');
ok($filter->execute, 'execute filter');
is(File::Compare::compare($expected_out_fastq, $out_fastq), 0, "fastq filtered as expected");

#print "$tmp_dir\n"; <STDIN>;
done_testing();
exit;

